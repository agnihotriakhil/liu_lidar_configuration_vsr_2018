
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <numeric>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Geometry>
using namespace std;
const double pi(3.14159265);

extern "C" 
struct coord{
    int x;
    int y;
    int z;
};

int manhattan_distance(coord a, coord b){
    return abs(a.x-b.x)+abs(a.y-b.y)+abs(a.z-b.z);
}

bool comp_xyz(const coord &a, const coord &b){
    if(a.x>b.x) return false;
    else if(a.x==b.x && a.y>b.y) return false;
    else if(a.x==b.x && a.y==b.y && a.z>b.z) return false;
    else return true;
}

bool comp_zxy(const coord &a, const coord &b){
    if(a.z>b.z) return false;
    else if(a.z==b.z && a.x>b.x) return false;
    else if(a.z==b.z && a.x==b.x && a.y>b.y) return false;
    else return true;
}

bool comp_yzx(const coord &a, const coord &b){
    if(a.y>b.y) return false;
    else if(a.y==b.y && a.z>b.z) return false;
    else if(a.y==b.y && a.z==b.z && a.x>b.x) return false;
    else return true;
}


void assign_cubes_to_subspaces(int* cube_num, float* cube_resolution, int* dead_zone_index, float* lidar_origin, int lidar_num,
                               int laser_num,  long subspace_num_max,   float* pitch_angle, vector<coord>* subspace_cube_coordinate){

    int cube_x_num=*(cube_num+0);
    int cube_y_num=*(cube_num+1);
    int cube_z_num=*(cube_num+2);
    float cube_resolution_x=*(cube_resolution+0);
    float cube_resolution_y=*(cube_resolution+1);
    float cube_resolution_z=*(cube_resolution+2);
    int dead_x_low=*(dead_zone_index+0);
    int dead_x_high=*(dead_zone_index+1);
    int dead_y_low=*(dead_zone_index+2);
    int dead_y_high=*(dead_zone_index+3);
    int dead_z_low=*(dead_zone_index+4);
    int dead_z_high=*(dead_zone_index+5);
  
    for(int i=0;i<cube_x_num;i++){
      for(int j=0;j<cube_y_num;j++){
        for(int k=0;k<cube_z_num;k++){
            //determine whether the cube is located in the dead zone
            if(i>=dead_x_low && i<=dead_x_high && j>=dead_y_low && j<=dead_y_high && k>=dead_z_low && k<=dead_z_high) continue;

            int subspace_index=0;
            for(int m=0;m<lidar_num;m++){
                int laser_index=0;
                float cube_x_world=cube_resolution_x*(i+0.5);
                float cube_y_world=cube_resolution_y*(j+0.5);
                float cube_z_world=cube_resolution_z*(k+0.5);
          
                // R P Y
                float roll=*(lidar_origin+6*m+3);
                float pitch=*(lidar_origin+6*m+4);
                float yaw=*(lidar_origin+6*m+5);
                Eigen::Vector3f euler_angle(yaw,pitch,roll);  
                Eigen::Matrix3f R;  
                R = Eigen::AngleAxisf(euler_angle[0], Eigen::Vector3f::UnitZ())  
                    * Eigen::AngleAxisf(euler_angle[1], Eigen::Vector3f::UnitY())  
                    * Eigen::AngleAxisf(euler_angle[2], Eigen::Vector3f::UnitX()); 

                Eigen::Vector3f v_3f;
                v_3f << cube_x_world,cube_y_world, cube_z_world;

                //transform from world coordinate to lidar local coordinate
                Eigen::Vector3f t;
                t<<*(lidar_origin+6*m),*(lidar_origin+6*m+1),*(lidar_origin+6*m+2);
                Eigen::Isometry3f T=Eigen::Isometry3f::Identity();  //Twl  lidar to world coordinate
                T.rotate(R);
                T.pretranslate(t);
                Eigen::Vector3f cube_local=T.inverse()*v_3f;        //Tlw  world to lidar coordinate
                
                float cube_x_local=cube_local[0];
                float cube_y_local=cube_local[1];
                float cube_z_local=cube_local[2];

                for(laser_index=0;laser_index<laser_num+1;laser_index++){
                    if(laser_index==laser_num)
                      break;
                    if(cube_z_local<=tan(pi*( *(pitch_angle+laser_index) )/180 )*sqrt(pow(cube_x_local,2)+pow(cube_y_local,2)) )
                      break;
                }
              
                subspace_index+=m*(laser_num+1)+laser_index;
            }

            // subspace_cube_num[subspace_index]+=1;
            coord tmp;
            tmp.x=i;
            tmp.y=j;
            tmp.z=k;
            subspace_cube_coordinate[subspace_index].push_back(tmp);
        }
      }
    }
}


extern "C" 
float subspace_segmentation(int* cube_num, float* cube_resolution, int* dead_zone_index, float* lidar_origin, int lidar_num,
                            int laser_num,  long subspace_num_max,   float* pitch_angle, float* result)
{
  
    // int* subspace_cube_num= new int[subspace_num_max]();  //={0};
    
    vector<coord>* subspace_cube_coordinate= new vector<coord>[subspace_num_max]();

    assign_cubes_to_subspaces(cube_num, cube_resolution, dead_zone_index, lidar_origin, lidar_num, laser_num, 
                            subspace_num_max, pitch_angle, subspace_cube_coordinate);
  
    vector<float> vsr; 
 
    for(int i=0;i<subspace_num_max;i++){

        if(subspace_cube_coordinate[i].size()==0)
            continue;
    
        vector<vector<coord>> results;
        int coord_num=subspace_cube_coordinate[i].size();
        int visited[coord_num]={0};
        int k=0;

        results.push_back(subspace_cube_coordinate[i]);
        // comment the above line if the LiDAR number is not 1
        
        for(int m=0;m<results.size();m++){
            vector<coord> tmp_cluster=results[m];
            float tmp_cube_num = float(tmp_cluster.size());

            //calculate surface_area and volume 
            int count=0;
            float surface_area=0;

            //xyz
            {
                sort(tmp_cluster.begin(),tmp_cluster.end(),comp_xyz);
                vector<coord>::iterator it = tmp_cluster.begin();
                coord cube_prev=*it;  //subspace_cube_coordinate[i].at(0);
                it++;
                for (; it != tmp_cluster.end(); it++){
                    if(it->x==cube_prev.x&&it->y==cube_prev.y) {
                      if(it->z==cube_prev.z+1)
                        count++;
                    }
                    cube_prev=*it;
                }
            }

            surface_area+=(2*tmp_cube_num-2*float(count))*cube_resolution_x*cube_resolution_y;
            count=0; 

            //zxy
            {
                sort(tmp_cluster.begin(),tmp_cluster.end(),comp_zxy);
                vector<coord>::iterator it = tmp_cluster.begin();
                coord cube_prev=*it;  //tmp_cluster.at(0);
                it++;
                for ( ;it != tmp_cluster.end(); it++){
                    if(it->z==cube_prev.z&&it->x==cube_prev.x){
                      if(it->y==cube_prev.y+1)
                        count++;
                    }
                    cube_prev=*it;
                }
            }

            surface_area+=(2*tmp_cube_num-2*float(count))*cube_resolution_z*cube_resolution_x;
            count=0; 

            //yzx
            {
                sort(tmp_cluster.begin(),tmp_cluster.end(),comp_yzx);
                vector<coord>::iterator it = tmp_cluster.begin();
                coord cube_prev=*it;  //tmp_cluster.at(0);
                it++;
                for ( ; it != tmp_cluster.end(); it++){
                    if(it->y==cube_prev.y&&it->z==cube_prev.z){
                      if(it->x==cube_prev.x+1)
                        count++;
                    }
                    cube_prev=*it;
                }
            }

            surface_area+=(2*tmp_cube_num-2*float(count))*cube_resolution_y*cube_resolution_z;
            count=0;

            //float surface_area=6*tmp_cube_num-2*count;
            float vsr_tmp=cube_resolution_x*cube_resolution_y*cube_resolution_z*float(tmp_cube_num)/surface_area;

            vsr.push_back(vsr_tmp);
        }
    }

    //find the maximum vsr
    float vsr_max=*max_element(vsr.begin(),vsr.end());
    float vsr_min=*min_element(vsr.begin(),vsr.end());
    
    float sum = std::accumulate(std::begin(vsr), std::end(vsr), 0.0);
    float mean = sum / vsr.size(); //均值

    float accum  = 0.0;
    
    std::for_each (std::begin(vsr), std::end(vsr), [&](const float d) {
      accum  += (d-mean)*(d-mean);
    });

    float stdev = sqrt(accum/(vsr.size()-1)); 

    delete[] subspace_cube_coordinate;
    delete[] subspace_cube_num;

    *(result+0)=vsr_max;
    *(result+1)=vsr_min;
    *(result+2)=stdev;
    *(result+3)=mean;

    return vsr_max;

}



