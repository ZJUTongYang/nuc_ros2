#pragma once
#include <vector>
#include <rclcpp/rclcpp.hpp>
#include <shape_msgs/msg/mesh.hpp>
#include <nuc_msgs/srv/get_nuc.hpp>
#include <nuc_msgs/srv/get_nuc_with_given_start.hpp>

#ifdef ENABLE_BENCHMARKING_3DCPP_INTERFACES
#include <benchmarking_3dcpp_interfaces/srv/get_nuc.hpp>
#endif

namespace nuc_ros2
{
    struct Facet
    {
    public: 
		// The index of the facet in the mesh
        int index_;

        Facet(){}
        Facet(const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver, 
            int tri_index, int subver_index, int child_in_parent_index)
        {
			child_.clear();
			index_ = tri_index;
			int v0 = mesh_tri[tri_index * 3];
			int v1 = mesh_tri[tri_index * 3 + 1];
			int v2 = mesh_tri[tri_index * 3 + 2];
			m_ = { (mesh_ver[v0 * 3] + mesh_ver[v1 * 3]) / 2, (mesh_ver[v0 * 3 + 1] + mesh_ver[v1 * 3 + 1]) / 2, (mesh_ver[v0 * 3 + 2] + mesh_ver[v1 * 3 + 2]) / 2 ,
				(mesh_ver[v1 * 3] + mesh_ver[v2 * 3]) / 2, (mesh_ver[v1 * 3 + 1] + mesh_ver[v2 * 3 + 1]) / 2, (mesh_ver[v1 * 3 + 2] + mesh_ver[v2 * 3 + 2]) / 2 ,
				(mesh_ver[v2 * 3] + mesh_ver[v0 * 3]) / 2, (mesh_ver[v2 * 3 + 1] + mesh_ver[v0 * 3 + 1]) / 2, (mesh_ver[v2 * 3 + 2] + mesh_ver[v0 * 3 + 2]) / 2
			};
			b_ = { (mesh_ver[v0 * 3] + mesh_ver[v1 * 3] + mesh_ver[v2 * 3]) / 3,
				(mesh_ver[v0 * 3 + 1] + mesh_ver[v1 * 3 + 1] + mesh_ver[v2 * 3 + 1]) / 3,
				(mesh_ver[v0 * 3 + 2] + mesh_ver[v1 * 3 + 2] + mesh_ver[v2 * 3 + 2]) / 3
			};
			polygon_center_ = { (mesh_ver[v1 * 3] + m_[1 * 3] + b_[0] + m_[0 * 3]) / 4, (mesh_ver[v1 * 3 + 1] + m_[1 * 3 + 1] + b_[0 + 1] + m_[0 * 3 + 1]) / 4, (mesh_ver[v1 * 3 + 2] + m_[1 * 3 + 2] + b_[0 + 2] + m_[0 * 3 + 2]) / 4,
								(mesh_ver[v2 * 3] + m_[2 * 3] + b_[0] + m_[1 * 3]) / 4, (mesh_ver[v2 * 3 + 1] + m_[2 * 3 + 1] + b_[0 + 1] + m_[1 * 3 + 1]) / 4, (mesh_ver[v2 * 3 + 2] + m_[2 * 3 + 2] + b_[0 + 2] + m_[1 * 3 + 2]) / 4,
								(mesh_ver[v0 * 3] + m_[0 * 3] + b_[0] + m_[2 * 3]) / 4, (mesh_ver[v0 * 3 + 1] + m_[0 * 3 + 1] + b_[0 + 1] + m_[2 * 3 + 1]) / 4, (mesh_ver[v0 * 3 + 2] + m_[0 * 3 + 2] + b_[0 + 2] + m_[2 * 3 + 2]) / 4
			};

			subver_index_ = subver_index;
			child_in_parent_index_ = child_in_parent_index;
        }
        void getCyclicCoveragePath(std::vector<int>& int_path, std::vector<double>& double_path) const
		{
			int subver_index = subver_index_ % 3;
			if (subver_index == 0)
			{
				int_path = { index_ * 3, index_ * 3 + 1, index_ * 3 + 2 };
				double_path.assign(polygon_center_.begin(), polygon_center_.end());
			}
			else if (subver_index == 1)
			{
				int_path = { index_ * 3 + 1, index_ * 3 + 2, index_ * 3 };
				double_path = { polygon_center_[3], polygon_center_[4], polygon_center_[5], polygon_center_[6], polygon_center_[7], polygon_center_[8], polygon_center_[0], polygon_center_[1], polygon_center_[2] };
			}
			else if (subver_index == 2)
			{
				int_path = { index_ * 3 + 2, index_ * 3, index_ * 3 + 1 };
				double_path = { polygon_center_[6], polygon_center_[7], polygon_center_[8], polygon_center_[0], polygon_center_[1], polygon_center_[2], polygon_center_[3], polygon_center_[4], polygon_center_[5] };
			}
		}

		std::vector<Facet*> child_;
		int child_in_parent_index_;
    private:

		std::vector<double> m_; // The Euclidean position of the midpoint of the three edges
		std::vector<double> b_; // The Euclidean position of the barocenter of the facet
		std::vector<double> polygon_center_; // The Euclidean position of the barocenter of sub-facets
		int subver_index_; // The order of the sub-facets designated by who the parent node is
    };

	class NUC: public rclcpp::Node
	{
		public:
			NUC();

			// void execute(const shape_msgs::msg::Mesh::ConstPtr& the_mesh);

			// bool executeService(nuc::GetNuc::Request& req, nuc::GetNuc::Response& resp);
			void getNUCCallback(const std::shared_ptr<nuc_msgs::srv::GetNuc::Request> request, 
				std::shared_ptr<nuc_msgs::srv::GetNuc::Response> response);

			void getNUCWithStartCallback(const std::shared_ptr<nuc_msgs::srv::GetNucWithGivenStart::Request> request, 
				std::shared_ptr<nuc_msgs::srv::GetNucWithGivenStart::Response> response);

		private:

			rclcpp::Service<nuc_msgs::srv::GetNucWithGivenStart>::SharedPtr nuc_with_start_service_;

			rclcpp::Service<nuc_msgs::srv::GetNuc>::SharedPtr nuc_service_;

			void convertMeshToVector(const shape_msgs::msg::Mesh& the_mesh, std::vector<int>& mesh_tri, std::vector<double>& mesh_ver);

#ifdef ENABLE_BENCHMARKING_3DCPP_INTERFACES
			rclcpp::Service<benchmarking_3dcpp_interfaces::srv::GetNuc>::SharedPtr nuc_benchmark_service_;
			void getNUCBenchmarkServiceCallback(const std::shared_ptr<benchmarking_3dcpp_interfaces::srv::GetNuc::Request> request, 
				std::shared_ptr<benchmarking_3dcpp_interfaces::srv::GetNuc::Response> response);
#endif
			
	};

};