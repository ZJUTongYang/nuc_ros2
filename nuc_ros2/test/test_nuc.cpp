#include <gtest/gtest.h>
#include <nuc_ros2/nuc.hpp>
#include <Eigen/Dense>

using namespace nuc_ros2;

TEST(NUC_Facet, CyclicCoveragePath_SingleTriangle)
{
    // Build a single triangle (3 vertices)
    std::vector<int> vertices_order = {0,1,2};
    std::vector<int> first_vertices_offset = {0, 3}; // one facet, offsets [0,3]

    std::vector<Eigen::VectorXd> vertex_positions;
    vertex_positions.resize(3);
    vertex_positions[0] = Eigen::VectorXd(3); vertex_positions[0] << 0.0, 0.0, 0.0;
    vertex_positions[1] = Eigen::VectorXd(3); vertex_positions[1] << 1.0, 0.0, 0.0;
    vertex_positions[2] = Eigen::VectorXd(3); vertex_positions[2] << 0.0, 1.0, 0.0;

    Facet f(vertices_order, vertex_positions, first_vertices_offset, 0);
    // initialize adjacency vector to indicate no neighbors
    f.adjacent_facets_index_.assign(f.num_vertices_, -1);
    // mark as root
    f.parent_index_ = ROOT_SENTINEL;
    f.parent_edge_index_ = 0;

    std::vector<Facet> all;
    all.push_back(f);

    auto result = all[0].getCyclicCoveragePath(all);
    auto int_path = result.first;
    auto double_path = result.second;

    // For a single triangle, we expect at least 3 subfacet centers
    EXPECT_GE(int_path.size(), 3u);
    EXPECT_EQ(int_path.size(), double_path.size());

    // Ensure the returned vectors contain valid coordinates
    for(const auto &v : double_path)
    {
        EXPECT_EQ(v.size(), 3);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
