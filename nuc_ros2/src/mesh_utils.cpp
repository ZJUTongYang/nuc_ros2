#include <nuc_ros2/mesh_utils.hpp>
#include <nuc_ros2/nuc.hpp>

namespace nuc_ros2
{

void convertMeshToVector(const shape_msgs::msg::Mesh& the_mesh, 
    std::vector<int>& facet_vertices, std::vector<double>& vertex_positions,
    unsigned int& coord_dim)
{
    unsigned int facet_num = the_mesh.triangles.size();
    unsigned int ver_num = the_mesh.vertices.size();
    // ROS Mesh vertices are 3D; we expose coord_dim so callers can support higher dimensions when available
    coord_dim = 3u;
    facet_vertices.resize(facet_num * 3);
    vertex_positions.resize(ver_num * coord_dim);

    #pragma omp parallel for
    for(size_t i = 0; i < facet_num; ++i)
    {
        facet_vertices[i*3] = the_mesh.triangles[i].vertex_indices[0];
        facet_vertices[i*3+1] = the_mesh.triangles[i].vertex_indices[1];
        facet_vertices[i*3+2] = the_mesh.triangles[i].vertex_indices[2];
    }

    #pragma omp parallel for
    for(size_t i = 0; i < ver_num; ++i)
    {
        vertex_positions[i*coord_dim + 0] = the_mesh.vertices[i].x;
        vertex_positions[i*coord_dim + 1] = the_mesh.vertices[i].y;
        vertex_positions[i*coord_dim + 2] = the_mesh.vertices[i].z;
    }
}

void compute_adjacency_matrix_directed(const std::vector<int>& vertices_order,
    const std::vector<int>& first_vertices_offset, 
    unsigned int ver_num,
    std::unordered_map<int, int>& amd)
{
    const unsigned int num_faces = first_vertices_offset.size();

    // 遍历网格中的每个面片。
    for (unsigned int face_idx = 0; face_idx < num_faces; ++face_idx)
    {

        // 获取当前面片在 vertices_order 中的起始索引。
        const int start_idx = first_vertices_offset[face_idx];

        // 确定当前面片在 vertices_order 中的结束索引。
        // 对于除最后一个面片之外的所有面片，其结束索引是下一个面片的起始索引。
        // 对于最后一个面片，其结束索引是 vertices_order 向量的末尾。
        int end_idx;
        if (face_idx < num_faces - 1) 
        {
            // exclusive end index
            end_idx = first_vertices_offset[face_idx + 1];
        }
        else 
        {
            end_idx = vertices_order.size();
        }

        // 计算当前面片的顶点数量。
        const int vertex_count = end_idx - start_idx;

        // 检查无效或退化面片。
        // 一个面片必须至少有3个顶点才能形成一个有效的多边形。
        if (vertex_count < 3)
        {
            printf("Degenerate face detected! This is not correct. \n");
            return ;
        }

        // 遍历当前多边形的边。
        for (int j = 0; j < vertex_count; ++j)
        {
            int u = vertices_order[start_idx + j];
            
            int v = vertices_order[start_idx + (j + 1) % vertex_count];

            int index = u + ver_num * v;

            amd[index] = face_idx;
        } // for each vertex
    }// for each facet
}

void determine_facet_adjacency(std::vector<nuc_ros2::Facet>& all_facets, const std::unordered_map<int, int>& amd, unsigned int ver_num)
{
    // This function determines the list of indices of adjacent facets for each facet
    for( auto& facet : all_facets)
    {
        std::vector<int>& vertices_indices = facet.vertices_indices_;
        facet.adjacent_facets_index_.clear();

        for (unsigned int i = 0; i < vertices_indices.size(); ++i)
        {
            int u = vertices_indices[i];
            int v = vertices_indices[(i + 1) % vertices_indices.size()];

            int reverse_edge_key = v + ver_num * u;
            if (amd.find(reverse_edge_key) != amd.end())
            {
                facet.adjacent_facets_index_.push_back(amd.at(reverse_edge_key));
            }
            else
            {
                facet.adjacent_facets_index_.push_back(-1);
            }
        }
    }
}

void assignFacetConnection(std::vector<Facet>& all_facets, int root_facet, int edge_index)
{
    // 深度优先扩展子节点，且保留邻接节点的循环索引顺序。
    all_facets[root_facet].parent_edge_index_ = edge_index;
    
    std::vector<int> adjacent_facet = all_facets[root_facet].adjacent_facets_index_;

    // edge_index指示的边不需要连回去，所以只有n-1条边需要连接
    for(unsigned int i = 1; i < adjacent_facet.size(); ++i)
    {
        int the_edge_to_connect = (edge_index + i) % adjacent_facet.size();
        int adj_index = adjacent_facet[the_edge_to_connect];
        // skip boundary edges
        if(adj_index == -1) continue;
        if(all_facets[adj_index].parent_index_ == NO_PARENT)
        {
            all_facets[adj_index].parent_index_ = static_cast<unsigned int>(root_facet);
            //我们需要确定对于child facet来说，其parent facet在哪个位置
            int loc = std::find(all_facets[adj_index].adjacent_facets_index_.begin(), 
                all_facets[adj_index].adjacent_facets_index_.end(), root_facet) - all_facets[adj_index].adjacent_facets_index_.begin();
            assignFacetConnection(all_facets, adj_index, loc);
        }
    }
    
}

} // namespace nuc_ros2