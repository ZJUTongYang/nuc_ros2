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
find_package(geometry_msgs REQUIRED)
find_package(tf2_eigen REQUIRED)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

# Optional dependency for benchmarking
find_package(benchmarking_3dcpp_interfaces)
if(benchmarking_3dcpp_interfaces_FOUND)
  message(STATUS "Found benchmarking_3dcpp_interfaces, enabling new message format")
  add_definitions(-DENABLE_BENCHMARKING_3DCPP_INTERFACES)
else()
  message(STATUS "benchmarking_3dcpp_interfaces not found, using old message format only")
endif()

include_directories(
  include 
  ${Eigen3_INCLUDE_DIRS}
  ${ament_INCLUDE_DIRS}
)


add_executable(nuc_node 
  src/nuc_node.cpp
  src/nuc.cpp)
target_include_directories(nuc_node PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
ament_target_dependencies(nuc_node PUBLIC rclcpp std_msgs shape_msgs nav_msgs geometry_msgs nuc_msgs tf2_eigen)

if(OpenMP_CXX_FOUND)
  target_link_libraries(nuc_node PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
endif()

if(benchmarking_3dcpp_interfaces_FOUND)
  ament_target_dependencies(nuc_node PUBLIC benchmarking_3dcpp_interfaces)
endif()

install(TARGETS
  nuc_node
  DESTINATION lib/${PROJECT_NAME})

add_library(nuc SHARED
  src/nuc.cpp)
target_include_directories(nuc PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
ament_target_dependencies(nuc PUBLIC
  rclcpp std_msgs shape_msgs nav_msgs geometry_msgs nuc_msgs tf2_eigen)

if(OpenMP_CXX_FOUND)
  target_link_libraries(nuc PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
endif()

if(benchmarking_3dcpp_interfaces_FOUND)
  ament_target_dependencies(nuc PUBLIC benchmarking_3dcpp_interfaces)
endif()


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
