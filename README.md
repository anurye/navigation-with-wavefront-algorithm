# Wavefront
## Introduction
The goal of this task is to introduce path planning implemented with the Wavefront Planner.

This algorithm is capable of finding the optimal route to the goal in the planning phase, in the initial state. For this purpose, it analyses the map and generates the path in an off-line manner. When finished, the robot simply follows the generated path.

In this task, you can run the simulation with the following command in the MATLAB console: 
`run_simulation(@solution6, false, [goal_x, goal_y], map_filename)`

- `solution6` is the control callback function you have to implement.
- `[goal_x, goal_y]` is the position of the goal point.
- `map_filename` is the name of the map file.

It also plots the state of the Wavefront Planner for debugging and visualization purposes.
The planner should run in the initial state. The generated collision-less trajectory should be executed using the controller implemented in the (Task 1). The robot should move to the subsequent cells on the path generated by Wavefront.  

## Using a map  
The planning algorithm requires a map of the environment. The map files are shipped together with the scene files `*.ttt`. Each scene file `<filename>.ttt` has its map file `<filename>.png`. You can use all scenes provided in the vrep_env directory to experiment and test your solution.    

To read a map file use the `imread` Matlab function, e.g.:  

`map = imread('vrep_env/map2.png')`

The map is a matrix, where each cell corresponds to the occupancy state of a particular place in the scene. The state of a cell is either 0 (occupied) or 255 (free), e.g. `map(1,1) == 0` means that the scene section with the lowest x and the lowest y is occupied. The map is indexed in ranges:

- for x: from 1 to size_x,
- for y: from 1 to size_y,

where `[size_y, size_x] = size(map);`. As the map corresponds to the scene, please take into account the units conversion, e.g. if the scene spans from -5m to 5m in x direction and from -10m to 10m in y direction, and the size of the map is 100 by 200 cells (`size(map)==[200,100]`), the units conversions are:
```
x_map = round(100 * ((x_world - (-5)) / (5 - (-5))))
y_map = round(200 * ((x_world - (-10)) / (10 - (-10))))
```
The inverse conversion can be easily calculated. Please note that the indexing of the map is (y,x).  
All environments span from -7.5m to 7.5m in both x and y directions, and the size of all maps is 100 by 100.  
Use plots to visualise the state of particle filter.   

# Setup
## Software components
Please read and follow the instruction carefully, to avoid errors and wasting your precious time!  

The programming/simulation environment used in EMOR tutorials reiles on four software components:  

- Matlab  
- Peter Corke’s RVC toolbox 
- The CoppeliaSim simulator (it was named V-REP in old times)  
- The Matlab bindings for CoppeliaSim
  
It is required to use a fairly recent version of Matlab. Please notice that Matlab versions older than 2011 may cause problems.  

## CoppeliaSim, Robotics Toolbox
The recommended procedure for setting up your computer is listed below:  
  
- Download the EMOR tutorials repository form github.com: the [ZIP snapshot](https://github.com/RCPRG-ros-pkg/emor_trs/archive/master.zip). Unzip the zip file 
  in the `~ws_emor` directory (in Ubuntu) or `C:/emor/ws_emor` directory (in Windows).  
- Download `CoppeliaSim EDU`. When you’re done downloading CoppeliaSim, unpack it in the `ws_emor` workspace directory.  
- Install the CoppeliaSim bindings for Matlab:  
1. You must copy three files from the directory of the CoppeliaSim app (downloaded at the previous step) to the directory named `youbot` within your local copy of the GitHub repository. The three files you need to copy are named  
- `remApi.m` – the file is located in `{CoppeliaSim_path}/programming/remoteApiBindings/matlab/matlab`; NOTE: in the newest version of CoppeliaSim, the files may be located in `{CoppeliaSim_path}/programming//legacyRemoteApi/remoteApiBindings/matlab/matlab`  
- `remoteApiProto.m` – the file is located in `{CoppeliaSim_path}/programming/remoteApiBindings/matlab/matlab`  
- `remoteApi.so` (if you use Linux) or `remoteApi.dll` (if you use Windows) or `remoteApi.dylib` (if you use a Mac). If you have a choice between a 32-bit or 64-bit `remoteApi`, pick the one that corresponds to your Matlab install (32-bit Matlab or 64-bit Matlab). If you have 32-bit Matlab, pick the 32-bit remoteApi, even if your kernel is 64-bit. The file is located in `{CoppeliaSim_path}/programming/remoteApiBindings/lib/lib` You will find these files in the directory containing the CoppeliaSim app. Look in the `programming/remoteApiBindings/lib/lib and programming/remoteApiBindings/matlab/matlab` subdirectories of the CoppeliaSim app directory (although this can change from version to version). You must copy these files to the youbot directory within your copy of the GitHub repo. If you closely followed the instructions above, the youbot directory is at `~/ws_emor/emor_trs/youbot` (Linux/Mac) or `C:/ws_emor/emor_trs/youbot` (Windows).
2. Run Matlab and change the current directory to the youbot directory (in Matlab Command Window), e.g. on Linux/Mac:

- `cd ~/ws_emor/emor_trs/youbot`
    
Then type (in Matlab):  
```
vrep=remApi('remoteApi');   
vrep.delete();
```  
If there is no error, the Matlab bindings to CoppeliaSim are loaded! If there is an error, check the steps above, and read CoppeliaSim Matlab bindings help.  
- In Matlab Command Window, run the `startup_robot.m` file provided via the Git repository. If you cloned the Git repository in the `ws_emor` workspace directory, run in Matlab Command Window:
 
`run('~/ws_emor/emor_trs/matlab/startup_robot.m');`

It will download and subsequently run the Peter Corke’s Robotics, Vision and Control toolkit. This line needs to be run every time Matlab restarts. You can add the script to your Matlab startup file.  


