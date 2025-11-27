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

namespace nuc_ros2
{
    void compute_adjacency_matrix_directed(const std::vector<int>& mesh_tri, unsigned int tri_num, int ver_num,
		std::unordered_map<int, int>& amd)
	{
		int index, first_ver, second_ver;
		for (unsigned int i = 0; i < tri_num; ++i)
		{
			if (mesh_tri[i * 3] == -1)
				continue;

			first_ver = mesh_tri[i * 3];
			second_ver = mesh_tri[i * 3 + 1];
			index = first_ver + ver_num * second_ver;
			amd[index] = i;

			first_ver = mesh_tri[i * 3 + 1];
			second_ver = mesh_tri[i * 3 + 2];
			index = first_ver + ver_num * second_ver;
			amd[index] = i;

			first_ver = mesh_tri[i * 3 + 2];
			second_ver = mesh_tri[i * 3];
			index = first_ver + ver_num * second_ver;
			amd[index] = i;
		}
	}

	Facet* insertTreeNode(Facet* parent, const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver, 
		int tri_index, int subver_index, int child_in_parent_index)
	{
		parent->child_.emplace_back(new Facet(mesh_tri, mesh_ver, tri_index, subver_index, child_in_parent_index));
		return parent->child_.back();
	}

	void deleteTreeBranch(Facet* node)
	{
		if (node == nullptr)
			return;
		for (Facet* child : node->child_)
		{
			deleteTreeBranch(child);
		}
		delete node;
	}

	void traverse(Facet* node, const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver, 
		std::vector<int>& topological_path, std::vector<double>& geometric_path, int position_to_insert)
	{
		if (node == nullptr)
			return;

		std::vector<int> cyclic_int_path;
		std::vector<double> cyclic_double_path;
		node->getCyclicCoveragePath(cyclic_int_path, cyclic_double_path);

		topological_path.insert(topological_path.begin() + position_to_insert, cyclic_int_path.begin(), cyclic_int_path.end());
		geometric_path.insert(geometric_path.begin() + position_to_insert * 3, cyclic_double_path.begin(), cyclic_double_path.end());

		int the_tri_index = node->index_;

		// int v0 = mesh_tri[the_tri_index * 3];
		// int v1 = mesh_tri[the_tri_index * 3+1];
		// int v2 = mesh_tri[the_tri_index * 3+2];
		int child_in_parent_index;
		int loc;
		for (auto iter = node->child_.begin(); iter != node->child_.end(); ++iter)
		{
			child_in_parent_index = (*iter)->child_in_parent_index_;
			if (child_in_parent_index == 0)
			{
				loc = std::find(topological_path.begin(), topological_path.end(), the_tri_index * 3 + 2) - topological_path.begin();
			}
			else if (child_in_parent_index == 1)
			{
				loc = std::find(topological_path.begin(), topological_path.end(), the_tri_index * 3) - topological_path.begin();
			}
			else if (child_in_parent_index == 2)
			{
				loc = std::find(topological_path.begin(), topological_path.end(), the_tri_index * 3 + 1) - topological_path.begin();
			}

			traverse(*iter, mesh_tri, mesh_ver, topological_path, geometric_path, loc + 1);
		}
	}

	std::pair<std::vector<int>, std::vector<double> > nuc(const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver, int initial_tri_index)
    {
		std::vector<int> topological_coverage_path;
		std::vector<double> geometric_coverage_path;

		Facet* root_ = nullptr;

		unsigned int tri_num = mesh_tri.size() / 3;
		unsigned int ver_num = mesh_ver.size() / 3;

		// We first collect the adjacency among facets
		std::unordered_map<int, int> amd; // <first_vertex_index + second_vertex_index * ver_num, facet_index>
		compute_adjacency_matrix_directed(mesh_tri, tri_num, ver_num, amd);

		// We avoid repetitive coverage
		std::vector<int> covered;
		covered.resize(tri_num, 0);

		// Invalid facets are marked as "covered", so that they are not considered during path deformation
		for (unsigned int i = 0; i < tri_num; ++i)
		{
			if (mesh_tri[i * 3] == -1)
			{
				covered[i] = 1;
			}
		}

		// The vertex indices of the source facet
		int v0 = mesh_tri[initial_tri_index * 3];
		int v1 = mesh_tri[initial_tri_index * 3 + 1];
		int v2 = mesh_tri[initial_tri_index * 3 + 2];

		// To avoid the problem in Theorem 1 in the paper, we select the best order of the sub-facets of the source facet
		// If any of the edges does not have adjacent facets, we select it as the "back"
		int subver_index;
		if (amd.find(v1 + v0 * ver_num) == amd.end())
		{
			subver_index = 0;
		}
		else if (amd.find(v2 + v1 * ver_num) == amd.end())
		{
			subver_index = 1;
		}
		else if (amd.find(v0 + v2 * ver_num) == amd.end())
		{
			subver_index = 2;
		}
		else
		{
			subver_index = 0;
		}
		
		// We put the initial facet at the root position
		root_ = new Facet(mesh_tri, mesh_ver, initial_tri_index, subver_index, -1);

		covered[initial_tri_index] = 1;
		std::list<Facet*> Q;
		Q.emplace_back(root_);

		int first_ver, second_ver;
		while (!Q.empty())
		{
			Facet* the_node = Q.front();
			Q.pop_front();

			int tri_index = the_node->index_;
			v0 = mesh_tri[tri_index * 3];
			v1 = mesh_tri[tri_index * 3 + 1];
			v2 = mesh_tri[tri_index * 3 + 2];

			first_ver = v0;
			second_ver = v1;
			std::unordered_map<int, int>::iterator iter = amd.find(second_ver + first_ver * ver_num);
			if (iter != amd.end())
			{
				int adj_tri_index = iter->second;
				if (covered[adj_tri_index] == 0)
				{
					if (mesh_tri[adj_tri_index * 3] == second_ver && mesh_tri[adj_tri_index * 3 + 1] == first_ver)
					{
						subver_index = 0;
					}
					else if (mesh_tri[adj_tri_index * 3 + 1] == second_ver && mesh_tri[adj_tri_index * 3 + 2] == first_ver)
					{
						subver_index = 1;
					}
					else if (mesh_tri[adj_tri_index * 3 + 2] == second_ver && mesh_tri[adj_tri_index * 3] == first_ver)
					{
						subver_index = 2;
					}
					else
					{
						subver_index = 0;
					}
					Facet* the_child = insertTreeNode(the_node, mesh_tri, mesh_ver, adj_tri_index, subver_index, 0);
					covered[adj_tri_index] = 1;
					Q.emplace_back(the_child);
				}
			}

			first_ver = v1;
			second_ver = v2;
			iter = amd.find(second_ver + first_ver * ver_num);
			if (iter != amd.end())
			{
				int adj_tri_index = iter->second;
				if (covered[adj_tri_index] == 0)
				{
					if (mesh_tri[adj_tri_index * 3] == second_ver && mesh_tri[adj_tri_index * 3 + 1] == first_ver)
					{
						subver_index = 0;
					}
					else if (mesh_tri[adj_tri_index * 3 + 1] == second_ver && mesh_tri[adj_tri_index * 3 + 2] == first_ver)
					{
						subver_index = 1;
					}
					else if (mesh_tri[adj_tri_index * 3 + 2] == second_ver && mesh_tri[adj_tri_index * 3] == first_ver)
					{
						subver_index = 2;
					}
					else
					{
						subver_index = 0;
					}
					Facet* the_child = insertTreeNode(the_node, mesh_tri, mesh_ver, adj_tri_index, subver_index, 1);
					covered[adj_tri_index] = 1;
					Q.emplace_back(the_child);
				}
			}

			first_ver = v2;
			second_ver = v0;
			iter = amd.find(second_ver + first_ver * ver_num);
			if (iter != amd.end())
			{
				int adj_tri_index = iter->second;
				if (covered[adj_tri_index] == 0)
				{
					if (mesh_tri[adj_tri_index * 3] == second_ver && mesh_tri[adj_tri_index * 3 + 1] == first_ver)
					{
						subver_index = 0;
					}
					else if (mesh_tri[adj_tri_index * 3 + 1] == second_ver && mesh_tri[adj_tri_index * 3 + 2] == first_ver)
					{
						subver_index = 1;
					}
					else if (mesh_tri[adj_tri_index * 3 + 2] == second_ver && mesh_tri[adj_tri_index * 3] == first_ver)
					{
						subver_index = 2;
					}
					else
					{
						subver_index = 0;
					}

					Facet* the_child = insertTreeNode(the_node, mesh_tri, mesh_ver, adj_tri_index, subver_index, 2);
					covered[adj_tri_index] = 1;
					Q.emplace_back(the_child);
				}
			}
		}

		// We report the geometric coverage path
		topological_coverage_path.clear();
		geometric_coverage_path.clear();

		int position_to_insert = 0;		
		traverse(root_, mesh_tri, mesh_ver, topological_coverage_path, geometric_coverage_path, position_to_insert);

		deleteTreeBranch(root_);

		return std::pair<std::vector<int>, std::vector<double> >(topological_coverage_path, geometric_coverage_path);
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

		return nuc(mesh_tri, mesh_ver, initial_tri_index);
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

		result = nuc(mesh_tri, mesh_ver, min_facet);

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

	void NUC::convertMeshToVector(const shape_msgs::msg::Mesh& the_mesh, 
		std::vector<int>& mesh_tri, std::vector<double>& mesh_ver)
	{
		unsigned int tri_num = the_mesh.triangles.size();
		unsigned int ver_num = the_mesh.vertices.size();
		mesh_tri.resize(tri_num*3);
		mesh_ver.resize(ver_num*3);

		#pragma omp parallel for
		for(size_t i = 0; i < tri_num; ++i)
		{
			mesh_tri[i*3] = the_mesh.triangles[i].vertex_indices[0];
			mesh_tri[i*3+1] = the_mesh.triangles[i].vertex_indices[1];
			mesh_tri[i*3+2] = the_mesh.triangles[i].vertex_indices[2];
		}

		#pragma omp parallel for
		for(size_t i = 0; i < ver_num; ++i)
		{
			mesh_ver[i*3] = the_mesh.vertices[i].x;
			mesh_ver[i*3+1] = the_mesh.vertices[i].y;
			mesh_ver[i*3+2] = the_mesh.vertices[i].z;
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
