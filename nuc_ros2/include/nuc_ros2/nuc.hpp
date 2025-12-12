#pragma once
#include <vector>
#include <rclcpp/rclcpp.hpp>
#include <shape_msgs/msg/mesh.hpp>
#include <geometry_msgs/msg/pose.hpp>
#include <nuc_msgs/srv/get_nuc.hpp>
#include <nuc_msgs/srv/get_nuc_with_given_start.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <tf2_eigen/tf2_eigen.hpp>

#ifdef ENABLE_BENCHMARKING_3DCPP_INTERFACES
#include <benchmarking_3dcpp_interfaces/srv/get_nuc.hpp>
#endif
namespace nuc_ros2
{
	struct Pose3D
	{
		// The homogeneous matrix
		Eigen::Matrix4d data;
	};

    struct Facet
    {
		// The basic structure for a facet in the given mesh 
		// TODO: add support for blended facets of different number of sides
		/*
	
		Vertices are in counter clockwise order
		m_i stays in the middle of v_i and v_{i+1}
		b is the barocenter

							  adj3
				v0-------------m3-------------v3
				| sub_facet 3  | sub_facet 2  |
		adj0	m0-------------b--------------m2   adj2
				| sub_facet 0  | sub_facet 1  |
				v1-------------m1-------------v2		
							  adj1
		
		*/
    public: 
		// The index of the facet in the mesh
        unsigned int index_;

        Facet(){}

		// This is for arbitrary blended mesh
		Facet(const std::vector<int>& vertices_order, const std::vector<Eigen::Vector3d>& vertices_position, 
			const std::vector<int>& first_vertices_offset, unsigned int facet_index)
		{
			// vertices order: [[v0, v1, v2, v3, v4, ...], [v0, v1, v2, v3, v4, ...], ...]
			// first_vertices_offset: [0, 3, 7, 11, ...]
			// vertices_position: [x0, y0, z0, x1, y1, z1, x2, y2, z2, ...
			// facet_index: the index of the facet in the mesh
			// subver_index: the index of the sub-facet in the facet. 
			// The index of a sub-facet is coincidentally identical to the index of its corresponding vertices in the vertices_order list (note that sub-facet i corresponds to the vertex i+1)
			// child_in_parent_index: the index of the child in the parent facet

			index_ = facet_index;
			parent_index_ = -1;

			int first_vertex_offset = first_vertices_offset[facet_index];
			int last_vertex_offset_exclusive;
			if(facet_index == first_vertices_offset.size() - 1)
			{
				last_vertex_offset_exclusive = vertices_order.size();
			}
			else
			{
				last_vertex_offset_exclusive = first_vertices_offset[facet_index + 1];
			}
			num_vertices_ = last_vertex_offset_exclusive - first_vertex_offset;

			b_ = Eigen::Vector3d(0, 0, 0);
			for(int i = first_vertex_offset; i < last_vertex_offset_exclusive; i++)
			{
				vertices_indices_.push_back(vertices_order[i]);
				v_.push_back(vertices_position[vertices_order[i]]);

				int vi = vertices_order[i];
				int viplus1;
				if (i == last_vertex_offset_exclusive - 1)
					viplus1 = vertices_order[first_vertex_offset];
				else
					viplus1 = vertices_order[i + 1];
				m_.push_back((vertices_position[vi] + vertices_position[viplus1]) / 2);
				b_ += m_.back();
			}
			b_ /= num_vertices_;

			polygon_center_.clear();
			for(size_t i = 0; i < num_vertices_; i++)
			{
				Eigen::Vector3d the_ver_pos;
				if (i == num_vertices_ - 1)
				{
					the_ver_pos = v_[0];
				}
				else
				{
					the_ver_pos = v_[i + 1];
				}

				Eigen::Vector3d sub_facet_center = (the_ver_pos + 
                                              m_[i] + 
                                              b_ + 
                                              m_[(i + 1) % num_vertices_]) / 4;
            
	            polygon_center_.push_back(sub_facet_center);
			}
			
			// This is the index of each sub-facet in the set of all sub-facets
			subfacet_indices_.resize(num_vertices_);
			for(size_t i = 0; i < num_vertices_; i++)
			{
				subfacet_indices_[i] = first_vertex_offset + i;
			}
		}

		// This function reports the sequence of sub-facets
        std::pair<std::vector<int>, std::vector<Eigen::Vector3d> > getCyclicCoveragePath(const std::vector<Facet>& all_facets) const
		{
			std::vector<int> int_path;
			std::vector<Eigen::Vector3d> double_path;
			
			// 循环路径是从parent_edge_index_这个subfacet开始，一直到parent_edge_index - 1结束
			for(unsigned int i = 0; i < num_vertices_; i++) 
			{
				int current_local_index = (parent_edge_index_ + i) % num_vertices_;

				int current_sub_facet_index = subfacet_indices_[current_local_index];
				int_path.push_back(current_sub_facet_index);
				double_path.push_back(polygon_center_[current_local_index]);
			}

			// 逐个检查adjacent node是否是child node，如果有child node，把它们的cyclic path加入到指定位置
			int position_to_insert = 1; // 如果接下来的第一个adjacent facet就是child facet且要返回路径，那么我们要把子路径插入到0元素和1元素之间
			for(size_t i = 1; i < num_vertices_; i++)
			{
				int current_local_index = (parent_edge_index_ + i) % num_vertices_;
				int adjacent_facet_index = adjacent_facets_index_[current_local_index];
				if(adjacent_facet_index != -1 && all_facets[adjacent_facet_index].parent_index_ == index_)
				{
					std::pair<std::vector<int>, std::vector<Eigen::Vector3d> > child_result = 
						all_facets[adjacent_facet_index].getCyclicCoveragePath(all_facets);
					int len = child_result.first.size();

					int_path.insert(int_path.begin()+position_to_insert, child_result.first.begin(), child_result.first.end());
					double_path.insert(double_path.begin()+position_to_insert, child_result.second.begin(), child_result.second.end());
					position_to_insert += len+1; // 跳过child path之后还得跳过一个自己的顶点
				}
				else
				{
					position_to_insert++;
				}
			}
			
			return std::make_pair(int_path, double_path);
		}

		int parent_index_; // -1 for the root facet
		int parent_edge_index_;

		unsigned int num_vertices_;
		std::vector<int> vertices_indices_;
		std::vector<int> subfacet_indices_;

		// m_i所在的边邻接的facet的索引
		std::vector<int> adjacent_facets_index_;
	
	private:

		std::vector<Eigen::Vector3d> v_; // The position of vertices

		std::vector<Eigen::Vector3d> m_; // The Euclidean position of the midpoint of the three edges

		Eigen::Vector3d b_; // The Euclidean position of the barocenter of the facet

		std::vector<Eigen::Vector3d> polygon_center_; // The Euclidean position of the barocenter of sub-facets

		int subver_index_; // The order of the sub-facets designated by who the parent node is

    };

	class NUC: public rclcpp::Node
	{
		public:
			NUC();

			void getNUCCallback(const std::shared_ptr<nuc_msgs::srv::GetNuc::Request> request, 
				std::shared_ptr<nuc_msgs::srv::GetNuc::Response> response);

			void getNUCWithStartCallback(const std::shared_ptr<nuc_msgs::srv::GetNucWithGivenStart::Request> request, 
				std::shared_ptr<nuc_msgs::srv::GetNucWithGivenStart::Response> response);

		private:

			rclcpp::Service<nuc_msgs::srv::GetNucWithGivenStart>::SharedPtr nuc_with_start_service_;

			rclcpp::Service<nuc_msgs::srv::GetNuc>::SharedPtr nuc_service_;

			void computeFacetNormals(const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver, 
				std::vector<double>& mesh_normal);

#ifdef ENABLE_BENCHMARKING_3DCPP_INTERFACES
			rclcpp::Service<benchmarking_3dcpp_interfaces::srv::GetNuc>::SharedPtr nuc_benchmark_service_;
			void getNUCBenchmarkServiceCallback(const std::shared_ptr<benchmarking_3dcpp_interfaces::srv::GetNuc::Request> request, 
				std::shared_ptr<benchmarking_3dcpp_interfaces::srv::GetNuc::Response> response);
#endif
			
	};
}
