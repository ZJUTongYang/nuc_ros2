cmake_minimum_required(VERSION 3.8)
project(nuc_ros2)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(shape_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(nuc_msgs REQUIRED)

include_directories(
  include 
  ${ament_INCLUDE_DIRS}
)


add_executable(nuc_node 
  src/nuc_node.cpp
  src/nuc.cpp)
target_include_directories(nuc_node PUBLIC
  ${ament_INCLUDE_DIRS}
)
ament_target_dependencies(nuc_node rclcpp std_msgs shape_msgs nav_msgs nuc_msgs)

install(TARGETS
  nuc_node
  DESTINATION lib/${PROJECT_NAME})

add_library(nuc SHARED
  src/nuc.cpp)
target_include_directories(nuc PUBLIC
  ${ament_INCLUDE_DIRS}
)
ament_target_dependencies(nuc PUBLIC
  rclcpp std_msgs shape_msgs nav_msgs nuc_msgs)

install(
  DIRECTORY include/
  DESTINATION include
)

install(
  TARGETS nuc
  EXPORT nuc
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
  )




if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
