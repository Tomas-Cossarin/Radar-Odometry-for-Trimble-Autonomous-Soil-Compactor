#ifndef COMPACTOR_H_INCLUDED
#define COMPACTOR_H_INCLUDED
#include "..\eigen-3.3.9\Eigen\Eigen"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

using namespace std;
using namespace Eigen;

class Radar
{
    public:
        Radar(double radar_pos_rax=0, double radar_pos_ray=0);

        double get_theta_ra(double depth, double theta_og);

        double x_pos_ra;
        double y_pos_ra;
};


class Epoch
{
    public:
        Epoch(int dataset, float ts=0, float heading=0, float x_velo=0, float y_velo=0);

        void get_odo_data(MatrixXd odo_csv, int ds);

        float ts;
        float heading;
        float x_velo;
        float y_velo;
        int dataset;
        VectorXd odo_data;
        float odo_x_pos;
        float odo_y_pos;
        float odo_heading;
        MatrixXd data;
};


class TestSuite
{
    public:
        TestSuite();

        void calibration_test(bool ransac=false);

        void ransac_test(bool create_outlier_data = false);

        void ls_only();

        bool prep_ransac;
        bool ransac;
        bool calibration;
        bool no_ransac;
        bool run_ls;
};


class Compactor
{
    public:
        Compactor(ofstream &logfile, int dataset, float x_pos=0, float y_pos=0, float x_pos_nh=0, float y_pos_nh=0,
                  float heading=0, float curr_ts=-1, float next_ts=-1, bool run_ls=false);

        void get_radar_data_col_lst();

        void new_pos(float speed, bool ls=false);

        void init_radar_positions();

        void log_epoch(Epoch epoch, ofstream &logfile);

        float x_pos;
        float y_pos;
        float heading;
        float x_pos_nh;
        float y_pos_nh;
        int dataset;
        TestSuite test_suite;
        float epoch_start;
        float curr_ts;
        float next_ts;
        float ls_heading;
        float ls_x_pos;
        float ls_y_pos;
        bool run_ls;
        vector<float> radar_orients;
        vector<Radar> radar_lst;
        vector<string> radar_data_col_lst;
};


#endif // COMPACTOR_H_INCLUDED
