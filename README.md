# Non-revisiting Uniform Coverage (NUC)

This repository is the ROS2 package for the journal paper entitled "Template-free Non-revisiting Uniform Coverage Path Planning on Curved Surfaces", IEEE/ASME Transactions on Mechatronics (T-Mech). 

[Paper Link](https://www.researchgate.net/publication/371391703)

This paper claimed that a physically uniform coverage path should be designed based on a physically uniform representation of the target surface, such as a triangular mesh, rather than enforcing a curvilinear coordinate. 

The supplementary code in this repo is capable of generating the non-repetitive coverage path on ARBITRARY connected triangular mesh. Input a uniform mesh, and you will get a uniform coverage. 

## Video Link

[Supplementary Video](https://drive.google.com/file/d/1sYnp-nKgyRzVhqUaI8ly20HRpIq9SC3B/view?usp=sharing)

## Dependencies

basic ROS2 installations

## Usage

1. Clone the code in the src folder of a ROS2 workspace
```
git clone https://github.com/ZJUTongYang/nuc_ros2.git
```

2. Compile the workspace
```
colcon build
```

This will create a ROS2 package, and a customized service package. 

3. Start the algorithm
```
ros2 run nuc_ros2 nuc_node
```
If a mesh is sent to the algorithm through service, a coverage path will be generated and returned. 

## Bug Report

The code does not have any randomness so bugs can be easily reproduced. If your surface triggers a bug in my code, please send the mesh data to me. 

## Cite This Paper
```
@article{Yang2023Template,
  title={Template-Free Nonrevisiting Uniform Coverage Path Planning on Curved Surfaces},
  author={Yang, Tong and Miro, Jaime Valls and Nguyen, Minh and Wang, Yue and Xiong, Rong},
  journal={IEEE/ASME Transactions on Mechatronics},
  year={2023},
  volume={28},
  number={4},
  pages={1853--1861},
  publisher={IEEE}
}
```


