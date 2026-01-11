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
		const std::vector<Eigen::VectorXd>& vertices_position, const std::vector<int>& first_vertices_offset, 
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
			unsigned int ver_num = vertices_position.size();
			compute_adjacency_matrix_directed(vertices_order, first_vertices_offset, ver_num, amd);

			determine_facet_adjacency(all_facets, amd, ver_num);

		int edge_index = find_best_starting_vertex(
			all_facets[initial_facet_index]);

	all_facets[initial_facet_index].parent_index_ = ROOT_SENTINEL; // so that the root facet won't be connected as a child
		assignFacetConnection(all_facets, initial_facet_index, edge_index);		

		// We report the geometric coverage path
		topological_coverage_path.clear();
		geometric_coverage_path.clear();

		std::pair<std::vector<int>, std::vector<Eigen::VectorXd> > child_result = 
			all_facets[initial_facet_index].getCyclicCoveragePath(all_facets);        
        
		topological_coverage_path = child_result.first;
		auto temp = child_result.second;
		for (auto &v : temp)
		{
			for(int d = 0; d < v.size(); ++d)
			{
				geometric_coverage_path.push_back(v[d]);
			}
		}

		return std::pair<std::vector<int>, std::vector<double> >(topological_coverage_path, geometric_coverage_path);							
 	}

	/**
	 * @brief 兼容性封装器，将旧的三角形网格接口适配到新的多边形网格接口。
	 *
	 * 此函数接受旧格式的网格数据，将其转换为新格式，调用新版本的 nuc_kernel，
	 * 并将结果转换回旧格式以保持兼容性。
	 *
	 * @param facet_vertices 旧格式的面片数据，每3个整数代表一个三角形（已用作兼容封装）。
	 * @param vertex_positions 旧格式的顶点数据，每3个浮点数代表一个顶点的坐标（已用作兼容封装）。
	 * @param initial_tri_index 初始面片的索引。
	 * @return std::pair 包含处理后的面片索引和顶点坐标（旧格式）。
	 */
	std::pair<std::vector<int>, std::vector<double> > nuc_kernel(const std::vector<int>& facet_vertices, 
		const std::vector<double>& vertex_positions, int initial_tri_index)
    {
		// TODO: 我们把它封装成任意混合多边形的形式
		std::vector<int> vertices_order;
		std::vector<int> first_vertices_offset;

	vertices_order.assign(facet_vertices.begin(), facet_vertices.end());

	const size_t num_facets = facet_vertices.size() / 3;
		first_vertices_offset.reserve(num_facets);
		for (size_t i = 0; i < num_facets; ++i) 
		{
			first_vertices_offset.push_back(i * 3);
		}

		// 转换顶点坐标: (std::vector<double>) -> (std::vector<Eigen::Vector3d>)
		const size_t num_vertices = vertex_positions.size() / 3;
		std::vector<Eigen::VectorXd> vertex_positions_eigen;
		vertex_positions_eigen.reserve(num_vertices);
		for (size_t i = 0; i < num_vertices; ++i) 
		{
			Eigen::VectorXd v(3);
			v << vertex_positions[i * 3], vertex_positions[i * 3 + 1], vertex_positions[i * 3 + 2];
			vertex_positions_eigen.push_back(v);
		}

		return nuc_kernel(vertices_order, vertex_positions_eigen, first_vertices_offset, initial_tri_index);
	}

	std::pair<std::vector<int>, std::vector<double> > nuc(const std::vector<int>& facet_vertices, const std::vector<double>& vertex_positions)
	{		
		// Here we allow for [-1, -1, -1] triangle facet, so we need to find the first valid facet
		int initial_tri_index = -1;
	unsigned int facet_num = facet_vertices.size()/3;
        
		for(unsigned int i = 0; i < facet_num; ++i)
		{
			if(facet_vertices[i*3] != -1)
			{
				initial_tri_index = i;
				break;
			}
		}

		if(initial_tri_index == -1)
			return std::pair<std::vector<int>, std::vector<double> >(std::vector<int>(), std::vector<double>());

	return nuc_kernel(facet_vertices, vertex_positions, initial_tri_index);
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
		
	std::vector<int> facet_vertices;
	std::vector<double> vertex_positions;
	convertMeshToVector(request->mesh, facet_vertices, vertex_positions);

	result = nuc(facet_vertices, vertex_positions);

		// The response of the benchmarking platform is a fully 3D one, so we need to assign the robot waypoint poses
		// The x-dir is the robot's heading direction, the y-dir is the robot's left direction, and the z-dir is the up direction
		// So for the case in the Yang2023Template paper, the tool's pointing direction is the -z axis
		std::vector<double> mesh_normal;
	unsigned int benchmark_coord_dim = 3;
	computeFacetNormals(facet_vertices, vertex_positions, mesh_normal, benchmark_coord_dim);

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
		
		std::vector<int> facet_vertices;
		std::vector<double> vertex_positions;
		convertMeshToVector(request->mesh, facet_vertices, vertex_positions);

		result = nuc(facet_vertices, vertex_positions);

		response->coverage.header.stamp = this->now();
		response->coverage.header.frame_id = request->frame_id;
		unsigned int num_waypoint = result.first.size();
		response->coverage.poses.resize(num_waypoint);

		// determine coordinate dimension from mesh vertices (compat: request->mesh has 3D vertices)
		unsigned int coord_dim = 3;
		if(request->mesh.vertices.size() > 0)
		{
			coord_dim = (vertex_positions.size() / request->mesh.vertices.size());
		}

		for(size_t i = 0; i < num_waypoint; ++i)
		{
			double x = 0, y = 0, z = 0;
			if(result.second.size() >= (i+1)*3)
			{
				x = result.second[i*coord_dim + 0];
				if(coord_dim > 1) y = result.second[i*coord_dim + 1];
				if(coord_dim > 2) z = result.second[i*coord_dim + 2];
			}
			response->coverage.poses[i].pose.position.x = x;
			response->coverage.poses[i].pose.position.y = y;
			response->coverage.poses[i].pose.position.z = z;
			response->coverage.poses[i].pose.orientation.w = 1;
		}
	}

	void NUC::getNUCWithStartCallback(const std::shared_ptr<nuc_msgs::srv::GetNucWithGivenStart::Request> request, 
		std::shared_ptr<nuc_msgs::srv::GetNucWithGivenStart::Response> response)
	{
		std::pair<std::vector<int>, std::vector<double> > result;
		
	std::vector<int> facet_vertices;
	std::vector<double> vertex_positions;
	convertMeshToVector(request->mesh, facet_vertices, vertex_positions);
		
		// We find the facet whose center is the closest to the start point
		double min_dist = 1000000;
		int min_facet = -1;
		unsigned int facet_num = facet_vertices.size()/3;
		unsigned int coord_dim = 3;
		if(request->mesh.vertices.size() > 0)
			coord_dim = vertex_positions.size() / request->mesh.vertices.size();

		for(size_t i = 0; i < facet_num; ++i)
		{
			// compute facet centroid in arbitrary dimension by averaging its vertex coordinates
			Eigen::VectorXd centroid = Eigen::VectorXd::Zero(coord_dim);
			int vi0 = facet_vertices[i*3 + 0];
			int vi1 = facet_vertices[i*3 + 1];
			int vi2 = facet_vertices[i*3 + 2];
			for(unsigned int d = 0; d < coord_dim; ++d)
			{
				double c0 = vertex_positions[vi0*coord_dim + d];
				double c1 = vertex_positions[vi1*coord_dim + d];
				double c2 = vertex_positions[vi2*coord_dim + d];
				centroid[d] = (c0 + c1 + c2) / 3.0;
			}

			double dx = request->start_pose.pose.position.x - (coord_dim>0 ? centroid[0] : 0.0);
			double dy = request->start_pose.pose.position.y - (coord_dim>1 ? centroid[1] : 0.0);
			double dz = request->start_pose.pose.position.z - (coord_dim>2 ? centroid[2] : 0.0);
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

	result = nuc_kernel(facet_vertices, vertex_positions, min_facet);

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


	void NUC::computeFacetNormals(const std::vector<int>& facet_vertices, const std::vector<double>& vertex_positions, 
				std::vector<double>& mesh_normal, unsigned int coord_dim)
	{
		// mesh_normal stores 3 components per facet (for compatibility with 3D consumers)
		mesh_normal.resize(facet_vertices.size());
		unsigned int facet_num = facet_vertices.size() / 3;

		#pragma omp parallel for
		for(size_t tri_index = 0; tri_index < facet_num; tri_index++)
		{
			// The vertices are listed in CCW order.
			int idx0 = facet_vertices[tri_index*3];
			int idx1 = facet_vertices[tri_index*3+1];
			int idx2 = facet_vertices[tri_index*3+2];

			Eigen::VectorXd v0 = Eigen::VectorXd::Zero(coord_dim);
			Eigen::VectorXd v1 = Eigen::VectorXd::Zero(coord_dim);
			Eigen::VectorXd v2 = Eigen::VectorXd::Zero(coord_dim);

			for(unsigned int d = 0; d < coord_dim; ++d)
			{
				v0[d] = vertex_positions[idx0*coord_dim + d];
				v1[d] = vertex_positions[idx1*coord_dim + d];
				v2[d] = vertex_positions[idx2*coord_dim + d];
			}

			Eigen::VectorXd edge1 = v1 - v0;
			Eigen::VectorXd edge2 = v2 - v0;

			Eigen::VectorXd normal_vec = Eigen::VectorXd::Zero(coord_dim);

			if(coord_dim == 3)
			{
				Eigen::Vector3d e1(edge1[0], edge1[1], edge1[2]);
				Eigen::Vector3d e2(edge2[0], edge2[1], edge2[2]);
				Eigen::Vector3d n = e1.cross(e2);
				double length = n.norm();
				if(length > 0) n /= length;
				normal_vec[0] = n[0]; normal_vec[1] = n[1]; normal_vec[2] = n[2];
			}
			else
			{
				// For higher dimensions, find a vector orthogonal to edge1 and edge2 via SVD
				Eigen::MatrixXd M(2, coord_dim);
				for(unsigned int d = 0; d < coord_dim; ++d)
				{
					M(0, d) = edge1[d];
					M(1, d) = edge2[d];
				}
				Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullV);
				Eigen::MatrixXd V = svd.matrixV();
				// take the last column of V (smallest singular value)
				normal_vec = V.col(coord_dim - 1);
				double normv = normal_vec.norm();
				if(normv > 0) normal_vec /= normv;
			}

			// Store first three components for compatibility with 3D consumers
			mesh_normal[tri_index*3+0] = (coord_dim > 0 ? normal_vec[0] : 0.0);
			mesh_normal[tri_index*3+1] = (coord_dim > 1 ? normal_vec[1] : 0.0);
			mesh_normal[tri_index*3+2] = (coord_dim > 2 ? normal_vec[2] : 0.0);
		}
	}

}
