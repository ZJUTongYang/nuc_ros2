#include <iostream>
#include <nuc_ros2/nuc.h>
#include <rclcpp/rclcpp.hpp>

int main(int argc, char** argv)
{
    rclcpp::init(argc, argv);

    auto node = std::make_shared<nuc_ros2::NUC>();

    rclcpp::spin(node);

    rclcpp::shutdown();

    return 0;

}