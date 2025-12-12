#include <vector>
#include <unordered_map>
#include <list>
#include <omp.h>
#include <nuc_ros2/nuc.hpp>
#include <iostream>
#include <algorithm>
#include <shape_msgs/msg/mesh.hpp>
#include <nav_msgs/msg/path.hpp>
#include <Eigen/Dense>
#include <nuc_ros2/mesh_utils.hpp>

namespace nuc_ros2
{

	/**
	 * @brief Determines the best starting vertex for a given face to avoid certain topological issues.
	 *
	 * This function searches for a boundary edge of the face. A boundary edge is an edge
	 * that is not shared with any other face. If a boundary edge is found, the function
	 * returns the index of the starting vertex of that edge. If no boundary edge is found,
	 * it returns 0 (the first vertex).
	 *
	 * @param initial_face_index The index of the face to analyze.
	 * @param vertices_order A flat vector of vertex indices for all faces.
	 * @param first_vertices_offset A vector of starting indices for each face in `vertices_order`.
	 * @param ver_num The total number of vertices in the mesh.
	 * @param amd The adjacency map, where keys are directed edges and values are face indices.
	 * @return The index of the best starting vertex within the face's vertex list.
	 */
	int find_best_starting_vertex(
		const Facet& the_facet)
	{
		std::vector<int> adjacent_facet = the_facet.adjacent_facets_index_;

		auto loc = std::find(adjacent_facet.begin(), adjacent_facet.end(), -1);

		if (loc != adjacent_facet.end())
		{
			// There is an edge that stays at the boundary of the surface
			return loc - adjacent_facet.begin();
		}
		else // this facet is in the middle of the surface
		{
			return 0;
		}
	}

	std::pair<std::vector<int>, std::vector<double> > nuc_kernel(const std::vector<int>& vertices_order, 
		const std::vector<Eigen::Vector3d>& vertices_position, const std::vector<int>& first_vertices_offset, 
		unsigned int initial_facet_index)
    {
		std::vector<int> topological_coverage_path;
		std::vector<double> geometric_coverage_path;
		
		unsigned int facet_num = first_vertices_offset.size();

		std::vector<Facet> all_facets;
		for(unsigned int i = 0; i < facet_num; ++i)
		{
			all_facets.emplace_back(Facet(vertices_order, vertices_position, first_vertices_offset, i));
		}
		// std::cout << "number of facets: " << all_facets.size() << std::endl;

		// We first collect the adjacency among facets
		std::unordered_map<int, int> amd; // <first_vertex_index + second_vertex_index * ver_num, facet_index>
		compute_adjacency_matrix_directed(vertices_order, first_vertices_offset, facet_num, amd);

		determine_facet_adjacency(all_facets, amd);

		int edge_index = find_best_starting_vertex(
			all_facets[initial_facet_index]);

		all_facets[initial_facet_index].parent_index_ = -2; // so that the root facet won't be connected as a child
		assignFacetConnection(all_facets, initial_facet_index, edge_index);		

		// We report the geometric coverage path
		topological_coverage_path.clear();
		geometric_coverage_path.clear();

		std::pair<std::vector<int>, std::vector<Eigen::Vector3d> > child_result = 
			all_facets[initial_facet_index].getCyclicCoveragePath(all_facets);		
		
		topological_coverage_path = child_result.first;
		auto temp = child_result.second;
		for (auto &v : temp)
		{
			geometric_coverage_path.push_back(v.x());
			geometric_coverage_path.push_back(v.y());
			geometric_coverage_path.push_back(v.z());
		}

		return std::pair<std::vector<int>, std::vector<double> >(topological_coverage_path, geometric_coverage_path);							
 	}

	/**
	 * @brief 兼容性封装器，将旧的三角形网格接口适配到新的多边形网格接口。
	 *
	 * 此函数接受旧格式的网格数据，将其转换为新格式，调用新版本的 nuc_kernel，
	 * 并将结果转换回旧格式以保持兼容性。
	 *
	 * @param mesh_tri 旧格式的面片数据，每3个整数代表一个三角形。
	 * @param mesh_ver 旧格式的顶点数据，每3个浮点数代表一个顶点的坐标。
	 * @param initial_tri_index 初始三角形的索引。
	 * @return std::pair 包含处理后的面片索引和顶点坐标（旧格式）。
	 */
	std::pair<std::vector<int>, std::vector<double> > nuc_kernel(const std::vector<int>& mesh_tri, 
		const std::vector<double>& mesh_ver, int initial_tri_index)
    {
		// TODO: 我们把它封装成任意混合多边形的形式
		std::vector<int> vertices_order;
		std::vector<int> first_vertices_offset;

		vertices_order.assign(mesh_tri.begin(), mesh_tri.end());

    	const size_t num_facets = mesh_tri.size() / 3;
		first_vertices_offset.reserve(num_facets);
		for (size_t i = 0; i < num_facets; ++i) 
		{
			first_vertices_offset.push_back(i * 3);
		}

		// 转换顶点坐标: (std::vector<double>) -> (std::vector<Eigen::Vector3d>)
    	const size_t num_vertices = mesh_ver.size() / 3;
		std::vector<Eigen::Vector3d> mesh_ver_eigen;
		mesh_ver_eigen.reserve(num_vertices);
		for (size_t i = 0; i < num_vertices; ++i) 
		{
			mesh_ver_eigen.emplace_back(
				mesh_ver[i * 3],    
				mesh_ver[i * 3 + 1],
				mesh_ver[i * 3 + 2] 
			);
		}

		return nuc_kernel(vertices_order, mesh_ver_eigen, first_vertices_offset, initial_tri_index);
	}

    std::pair<std::vector<int>, std::vector<double> > nuc(const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver)
	{		
		// Here we allow for [-1, -1, -1] triangle facet, so we need to find the first valid facet
		int initial_tri_index = -1;
		unsigned int tri_num = mesh_tri.size()/3;
		
		for(unsigned int i = 0; i < tri_num; ++i)
		{
			if(mesh_tri[i*3] != -1)
			{
				initial_tri_index = i;
				break;
			}
		}

		if(initial_tri_index == -1)
			return std::pair<std::vector<int>, std::vector<double> >(std::vector<int>(), std::vector<double>());

		return nuc_kernel(mesh_tri, mesh_ver, initial_tri_index);
	}

	NUC::NUC(): Node("nuc")
	{

		nuc_service_ = this->create_service<nuc_msgs::srv::GetNuc>("get_nuc", 
			std::bind(&NUC::getNUCCallback, this, std::placeholders::_1, std::placeholders::_2));
	
		nuc_with_start_service_ = this->create_service<nuc_msgs::srv::GetNucWithGivenStart>("get_nuc_with_start", 
			std::bind(&NUC::getNUCWithStartCallback, this, std::placeholders::_1, std::placeholders::_2));

#ifdef ENABLE_BENCHMARKING_3DCPP_INTERFACES
			nuc_benchmark_service_ = this->create_service<benchmarking_3dcpp_interfaces::srv::GetNuc>("get_Yang2023Template_benchmark", 
				std::bind(&NUC::getNUCBenchmarkServiceCallback, this, std::placeholders::_1, std::placeholders::_2));
#endif
	}

#ifdef ENABLE_BENCHMARKING_3DCPP_INTERFACES
	void NUC::getNUCBenchmarkServiceCallback(const std::shared_ptr<benchmarking_3dcpp_interfaces::srv::GetNuc::Request> request,
		std::shared_ptr<benchmarking_3dcpp_interfaces::srv::GetNuc::Response> response)
	{
	    std::pair<std::vector<int>, std::vector<double> > result;
		
		std::vector<int> mesh_tri;
		std::vector<double> mesh_ver;
		convertMeshToVector(request->mesh, mesh_tri, mesh_ver);

		result = nuc(mesh_tri, mesh_ver);

		// The response of the benchmarking platform is a fully 3D one, so we need to assign the robot waypoint poses
		// The x-dir is the robot's heading direction, the y-dir is the robot's left direction, and the z-dir is the up direction
		// So for the case in the Yang2023Template paper, the tool's pointing direction is the -z axis
		std::vector<double> mesh_normal;
		computeFacetNormals(mesh_tri, mesh_ver, mesh_normal);

		// By having the subfacets index, we know the facet that the waypoint stays in
		std::vector<int> subfacet_index_path = result.first;
		std::vector<double> waypoint_path = result.second;

		std::vector<Pose3D> path_3d;
		path_3d.resize(subfacet_index_path.size());

		// We set up the xyz
		// The waypoint is set as a [16] array
		// position is copied as is
		// z-direction (outer normal of the facet)
		#pragma omp parallel for
		for(size_t i = 0; i < path_3d.size(); ++i)
		{
			auto& data = path_3d[i].data;

			data(0, 3) = waypoint_path[i*3];
			data(1, 3) = waypoint_path[i*3+1];
			data(2, 3) = waypoint_path[i*3+2];
			data(3, 3) = 1;

			int facet_index = subfacet_index_path[0] / 3; 

			data(0, 2) = mesh_normal[facet_index*3];
			data(1, 2) = mesh_normal[facet_index*3+1];
			data(2, 2) = mesh_normal[facet_index*3+2];			
		}


		// We estimate the X-direction (the diff of the beginning two waypoint)
		for(size_t i = 0; i < path_3d.size(); ++i)
		{
			Eigen::Vector3d curr, prev, next;
			curr << waypoint_path[i*3], waypoint_path[i*3+1], waypoint_path[i*3+2];
			if(i == 0)
			{
				prev = curr;
			}
			else
			{
				prev << waypoint_path[(i-1)*3], waypoint_path[(i-1)*3+1], waypoint_path[(i-1)*3+2];
			}

			if(i == path_3d.size()-1)
			{
				next = curr;
			}
			else
			{
				next << waypoint_path[(i+1)*3], waypoint_path[(i+1)*3+1], waypoint_path[(i+1)*3+2];
			}

			Eigen::Vector3d mean_diff = (next - prev) / 2.0;

			auto& data = path_3d[i].data;
			data(0, 0) = mean_diff(0);
			data(1, 0) = mean_diff(1);
			data(2, 0) = mean_diff(2);

			// We set the y-dir as the cross product of x-dir and z-dir
			data(0, 1) = data(1, 2) * data(2, 0) - data(2, 2) * data(1, 0);
			data(1, 1) = data(2, 2) * data(0, 0) - data(0, 2) * data(2, 0);
			data(2, 1) = data(0, 2) * data(1, 0) - data(1, 2) * data(0, 0);

			// We normalize each direction
			data.col(0).normalize();
			data.col(1).normalize();
			data.col(2).normalize();
		}

		response->coverage.poses.resize(path_3d.size());
		for(size_t i = 0; i < path_3d.size(); ++i)
		{
			response->coverage.poses[i].pose.position.x = path_3d[i].data(0, 3);
			response->coverage.poses[i].pose.position.y = path_3d[i].data(1, 3);
			response->coverage.poses[i].pose.position.z = path_3d[i].data(2, 3);

			// We transform the rotation matrix to quaternion
			Eigen::Quaterniond q(path_3d[i].data.block<3, 3>(0, 0));
			response->coverage.poses[i].pose.orientation.x = q.x();
			response->coverage.poses[i].pose.orientation.y = q.y();
			response->coverage.poses[i].pose.orientation.z = q.z();
			response->coverage.poses[i].pose.orientation.w = q.w();
		}

	}
#endif

	void NUC::getNUCCallback(const std::shared_ptr<nuc_msgs::srv::GetNuc::Request> request, 
		std::shared_ptr<nuc_msgs::srv::GetNuc::Response> response)
	{
		std::pair<std::vector<int>, std::vector<double> > result;
		
		std::vector<int> mesh_tri;
		std::vector<double> mesh_ver;
		convertMeshToVector(request->mesh, mesh_tri, mesh_ver);

		result = nuc(mesh_tri, mesh_ver);

		response->coverage.header.stamp = this->now();
		response->coverage.header.frame_id = request->frame_id;
		unsigned int num_waypoint = result.first.size();
		response->coverage.poses.resize(num_waypoint);

		for(size_t i = 0; i < num_waypoint; ++i)
		{
			response->coverage.poses[i].pose.position.x = result.second[i*3];
			response->coverage.poses[i].pose.position.y = result.second[i*3+1];
			response->coverage.poses[i].pose.position.z = result.second[i*3+2];
			response->coverage.poses[i].pose.orientation.w = 1;
		}
	}

	void NUC::getNUCWithStartCallback(const std::shared_ptr<nuc_msgs::srv::GetNucWithGivenStart::Request> request, 
		std::shared_ptr<nuc_msgs::srv::GetNucWithGivenStart::Response> response)
	{
		std::pair<std::vector<int>, std::vector<double> > result;
		
		std::vector<int> mesh_tri;
		std::vector<double> mesh_ver;
		convertMeshToVector(request->mesh, mesh_tri, mesh_ver);
		
		// We find the facet whose center is the closest to the start point
		double min_dist = 1000000;
		int min_facet = -1;
		unsigned int tri_num = mesh_tri.size()/3;
		for(size_t i = 0; i < tri_num; ++i)
		{
			double dx = request->start_pose.pose.position.x - mesh_ver[i*3];
			double dy = request->start_pose.pose.position.y - mesh_ver[i*3+1];
			double dz = request->start_pose.pose.position.z - mesh_ver[i*3+2];
			double dist = dx*dx + dy*dy + dz*dz;
			if(dist < min_dist)
			{
				min_dist = dist;
				min_facet = i;
			}
		}
		if(min_facet == -1)
		{
			std::cout << "nuc: cannot find the closest facet" << std::endl;
			return;
		}

		result = nuc_kernel(mesh_tri, mesh_ver, min_facet);

		response->coverage.header.stamp = this->now();
		response->coverage.header.frame_id = request->frame_id;
		unsigned int num_waypoint = result.first.size();
		response->coverage.poses.resize(num_waypoint);

		for(size_t i = 0; i < num_waypoint; ++i)
		{
			response->coverage.poses[i].pose.position.x = result.second[i*3];
			response->coverage.poses[i].pose.position.y = result.second[i*3+1];
			response->coverage.poses[i].pose.position.z = result.second[i*3+2];
			response->coverage.poses[i].pose.orientation.w = 1;
		}
	}


	void NUC::computeFacetNormals(const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver, 
				std::vector<double>& mesh_normal)
	{
		mesh_normal.resize(mesh_tri.size());
		unsigned int tri_num = mesh_tri.size() / 3;

		#pragma omp parallel for
		for(size_t tri_index = 0; tri_index < tri_num; tri_index++)
		{
			// The vertices are listed in CCW order.
			int idx0 = mesh_tri[tri_index*3];
			int idx1 = mesh_tri[tri_index*3+1];
			int idx2 = mesh_tri[tri_index*3+2];

			double v0[3] = {mesh_ver[idx0*3], mesh_ver[idx0*3+1], mesh_ver[idx0*3+2]};
			double v1[3] = {mesh_ver[idx1*3], mesh_ver[idx1*3+1], mesh_ver[idx1*3+2]};
			double v2[3] = {mesh_ver[idx2*3], mesh_ver[idx2*3+1], mesh_ver[idx2*3+2]};

			double edge1[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
	        double edge2[3] = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};
        
			// e1 cross_product e2 is the outer normal vector
			double normal[3];
			normal[0] = edge1[1]*edge2[2] - edge1[2]*edge2[1];
			normal[1] = edge1[2]*edge2[0] - edge1[0]*edge2[2];
			normal[2] = edge1[0]*edge2[1] - edge1[1]*edge2[0];

			double length = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
			if(length > 0) {
				normal[0] /= length;
				normal[1] /= length;
				normal[2] /= length;
			}
			else
			{
				std::cout << "Error: a facet has normal vector of length 0" << std::endl;
			}

			mesh_normal[tri_index*3+0] = normal[0];
			mesh_normal[tri_index*3+1] = normal[1];
			mesh_normal[tri_index*3+2] = normal[2];
		}
	}

}
