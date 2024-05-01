#include <iostream>
#include "..\eigen-3.3.9\Eigen\Eigen"
#include <fstream>
#include <string>
#include <vector>

#include "classes.h"
#include "functions.h"

using namespace std;
using namespace Eigen;

void run_odo(Compactor compactor, ofstream& logfile);
void main_loop(Compactor compactor, MatrixXd radar_csv, ofstream& logfile);
void plot_odo_results(Compactor compactor, vector<float> cal_adjusts);


void run_odo(Compactor compactor, ofstream& logfile)
{
    //vector<float> cal_adjusts = {0,0,0,0,0,0};
    vector<float> cal_adjusts = {-0.1428571428571429, 0, 0, -0.020408163265306145, -0.1428571428571429, 0};
    //cal_adjusts = {-0.08333333333333337, 0.08333333333333326, 0, -0.9166666666666666, 0, 0};
    //cal_adjusts = {1.0, 0, 0.0, 0, -0.16666666666666674, 0};
    //cal_adjusts = {0.33333333333333326, 0, 0, 0, -0.8333333333333334, 0};

    compactor.radar_orients = get_radar_orient(cal_adjusts);

    string radar_data_path  = "../data/rad_data_ds_c" + to_string(compactor.dataset) + ".csv";
    MatrixXd radar_csv = read_csv(radar_data_path);

    main_loop(compactor, radar_csv, logfile);
    plot_odo_results(compactor, cal_adjusts);
}


void main_loop(Compactor compactor, MatrixXd radar_csv, ofstream& logfile)
{
    cout<<"MAIN LOOP"<<endl;
    string odo_data_path = "../data/odo_data_ds" + to_string(compactor.dataset) + ".csv";
    MatrixXd odo_csv = read_csv(odo_data_path);

    float mid_ts;
    int counter = 0;
    int cols = radar_csv.cols();

    for(int i=0; i<radar_csv.rows(); i++)
    {
        VectorXd line = radar_csv.row(i);
        compactor.next_ts = line[1];

        if(compactor.curr_ts == -1)
        {
            compactor.curr_ts = line[1];
            mid_ts = compactor.curr_ts;
        }
        if(compactor.next_ts != mid_ts)
        {
            MatrixXd epoch_data = radar_csv.block(i-counter, 0, counter, cols);
            vector<float> nh = process_epoch(epoch_data, compactor, odo_csv, logfile);

            compactor.x_pos_nh += nh[0];
            compactor.y_pos_nh += nh[1];

            compactor.curr_ts = compactor.next_ts;
            mid_ts = compactor.next_ts;
            counter = 1;
        }
        else
            counter++;
    }
    logfile.close();
}


void plot_odo_results(Compactor compactor, vector<float> cal_adjusts)
{
    string ground_truth_path = "../data/consolidated_data_gt-ds-de" + to_string(compactor.dataset) + ".csv";
    MatrixXd ground_truth_csv = read_csv(ground_truth_path);

    string curr_data_path = "../logfiles/consolidated_data_ds" + to_string(compactor.dataset) + ".csv";
    MatrixXd curr_data_csv = read_csv(curr_data_path);

    VectorXd test_x_nh = curr_data_csv.col(7);
    VectorXd test_y_nh = curr_data_csv.col(8);

    VectorXd odo_x0 = ground_truth_csv.col(11);
    VectorXd odo_y0 = ground_truth_csv.col(12);
    VectorXd odo_x = VectorXd::Zero(test_x_nh.size());
    VectorXd odo_y = VectorXd::Zero(test_y_nh.size());
    int first_x = odo_x0(0);
    int first_y = odo_y0(0);

    for(int r=0; r<test_x_nh.size(); r++)
    {
        odo_x(r) = odo_x0(r) - first_x;
        odo_y(r) = odo_y0(r) - first_y;
    }

    vector<float> pos_neg_area_x = find_error(odo_x, test_x_nh);
    cout<<"Positive area between calculated x and GNSS x: "<<pos_neg_area_x[0]<<endl;
    cout<<"Negative area between calculated x and GNSS x: "<<pos_neg_area_x[1]<<endl;
    cout<<"Total x error: "<<abs(pos_neg_area_x[0])+abs(pos_neg_area_x[1])<<endl;

    vector<float> pos_neg_area_y = find_error(odo_y, test_y_nh);
    cout<<"Positive area between calculated y and GNSS y: "<<pos_neg_area_y[0]<<endl;
    cout<<"Negative area between calculated y and GNSS y: "<<pos_neg_area_y[1]<<endl;
    cout<<"Total y error: "<<abs(pos_neg_area_y[0])+abs(pos_neg_area_y[1])<<endl;

    VectorXd time = VectorXd::Zero(ground_truth_csv.rows());
    float start_time = (int) ground_truth_csv(0,0);

    for(int i=0; i<ground_truth_csv.rows(); i++)
        time(i) = ground_truth_csv(i,0) - start_time;

    bool calibrated = false;
    for(int c=0; c<6; c++)
    {
        if(cal_adjusts[c] != 0)
        {
            calibrated = true;
            cout<<"Calibrated"<<endl;
            break;
        }
    }
    if(calibrated==false)
        cout<<endl<<"Uncalibrated"<<endl;

    MatrixXd dep_varsx = MatrixXd::Zero(odo_x.size(), 2);
    MatrixXd dep_varsy = MatrixXd::Zero(odo_x.size(), 2);
    dep_varsx <<odo_x, test_x_nh;
    dep_varsy <<odo_y, test_y_nh;

    MatrixXd dep_vars = MatrixXd::Zero(odo_x.size(), 4);
    dep_vars <<dep_varsx, dep_varsy;

    MatrixXd time_dep_vars = MatrixXd::Zero(odo_x.size(), 5);
    time_dep_vars<<time, dep_vars;

    ofstream out("../Output_ds"+ to_string(compactor.dataset) +".txt");
    out<<time_dep_vars;
    out.close();
}


int main()
{
    int dataset = 2;

    cout<<"START"<<endl<<endl;

    string logfile_path = "../logfiles/consolidated_data_ds"+ to_string(dataset) +".csv";
    ofstream logfile;
    logfile.open(logfile_path);

    Compactor compactor (logfile, dataset);
    compactor.init_radar_positions();

    compactor.test_suite.ls_only();   //uses LS without ransac

    run_odo(compactor, logfile);

    return 0;
}
