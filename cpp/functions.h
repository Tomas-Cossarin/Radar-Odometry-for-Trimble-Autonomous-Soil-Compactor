#ifndef PROCESS_EPOCH_H_INCLUDED
#define PROCESS_EPOCH_H_INCLUDED
#include "..\eigen-3.3.9\Eigen\Eigen"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "classes.h"

using namespace std;
using namespace Eigen;

void write_to_log(ofstream &logfile, string str);

MatrixXd read_csv(string &filepath);

vector<float> find_error(VectorXd odo, VectorXd test_nh);

float trapz(VectorXd y);

//vector<vector<vector<float>>> get_odo_data_ts_dict(MatrixXd odo_csv, int ds=1);

VectorXd get_odo_ts_data(float ts, MatrixXd odo_csv);

float get_odo_heading_ts(VectorXd odo_ts_data, int ds);

vector<float> quat_to_ang(float x, float y, float z, float w);

vector<float> process_epoch(MatrixXd epoch_data, Compactor compactor, MatrixXd odo_csv, ofstream &logfile);

vector<float> new_pos_nh(Compactor compactor, float v_x, float v_y);

vector<float> get_radar_orient(vector<float> cal_adjusts={0,0,0,0,0,0}, int radar_idx=-1, float adjust=0);

vector<float> get_epoch_velo(Compactor compactor, MatrixXd epoch_data, bool ls=false);

vector<vector<vector<double>>> get_radar_dict(MatrixXd epoch_data);

MatrixXd create_matrices(Compactor compactor, vector<vector<vector<double>>> radar_dict);

vector<float> least_squares_naive(MatrixXd A, VectorXd speed_radials);

vector<float> least_squares_ransac(MatrixXd A, VectorXd speed_radials);


#endif // PROCESS_EPOCH_H_INCLUDED
