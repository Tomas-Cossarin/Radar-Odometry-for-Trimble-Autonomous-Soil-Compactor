#include "functions.h"
#include "classes.h"

void write_to_log(ofstream &logfile, string str)
{
    logfile << str;
    logfile << '\n';
}


MatrixXd read_csv(string &filepath)
{
    ifstream file0(filepath);
    float numRows = count(istreambuf_iterator<char>(file0), istreambuf_iterator<char>(), '\n');
    file0.close();

    ifstream file1(filepath);
    string cell1;
    string line1;
    getline(file1, line1);
    stringstream lineStream1(line1);
    int numCols = 0;
    while(getline(lineStream1, cell1, ','))
        numCols++;
    file1.close();

    MatrixXd out = MatrixXd::Zero(numRows-1, numCols);
    ifstream file(filepath);

    for(int i=0; i<numRows; i++)
    {
        vector<string> row;
        string cell;
        string line;
        getline(file, line);
        stringstream lineStream(line);

        while(getline(lineStream, cell, ','))
        {
            row.push_back(cell);
        }
        if(i!=0)
        {
            for(int j=0; j<row.size(); j++)
            {
                out(i-1, j) = stof(row[j]);
            }
        }
    }
    file.close();
    return out;
}


vector<float> find_error(VectorXd odo, VectorXd test_nh)
{
    VectorXd diff = odo - test_nh;
    VectorXd posPart = VectorXd::Zero(diff.size());
    VectorXd negPart = VectorXd::Zero(diff.size());
    for(int i=0; i<diff.size(); i++)
    {
        float di = diff(i);
        if(di > 0)
            posPart(i) = di;

        else if(di < 0)
            negPart(i) = -di;
    }

    vector<float> pos_neg_area = {trapz(posPart), trapz(negPart)};
    return pos_neg_area;
}


float trapz(VectorXd y)
{
    float area = 0;
    for(int i=0; i<y.size()-1; i++)
    {
        area += (y[i+1] + y[i]);
    }
    return area/float(2);
}


/*vector<vector<vector<float>>> get_odo_data_ts_dict(MatrixXd odo_csv, int ds)
{
    int rows = (int)odo_csv(odo_csv.rows()-1,0) - (int)odo_csv(0,0) + 1;
    int aisles = odo_csv.cols();
    vector<vector<vector<float>>> odo_data_ts_dict;
    for(int j=0; j<rows; j++)
        odo_data_ts_dict.push_back({});

    int prev_time = 0;
    int j = -1;
    for(int i=0; i<odo_csv.rows(); i++)
    {
        vector<float> col;
        for(int k=0; k<aisles; k++)
            col.push_back(odo_csv.row(i)[k]);

        int time = (int) odo_csv(i,0);
        if(time != prev_time)
            j++;
        prev_time = time;

        int n = odo_data_ts_dict[j].size();
        odo_data_ts_dict[j].push_back({});
        for(int k=0; k<aisles; k++)
            odo_data_ts_dict[j][n].push_back(col[k]);
    }
    return odo_data_ts_dict;
}*/


VectorXd get_odo_ts_data(float ts, MatrixXd odo_csv)
{
    for(int i=0; i<odo_csv.rows()-1; i++)
    {
        float tsi = odo_csv(i,0);
        if(tsi == ts)
            return odo_csv.row(i);

        else if(tsi > ts)
        {
            if( abs(tsi-odo_csv(i-1,0)) < abs(tsi-odo_csv(i,0)) )
                return odo_csv.row(i-1);
            else
                return odo_csv.row(i);
        }
    }
}


float get_odo_heading_ts(VectorXd odo_ts_data, int ds)
{
    float qx = odo_ts_data[4];
    float qy = odo_ts_data[5];
    float qz = odo_ts_data[6];
    float qw = odo_ts_data[7];

    vector<float> xyz = quat_to_ang(qx, qy, qz, qw);
    return xyz[2]*M_PI/180;
}


vector<float> quat_to_ang(float x, float y, float z, float w)
{
    // Params: quaternion
    // Returns: x, y, z in degrees

    x = atan2((2.0*(w*x + y*z)), (1.0 - 2.0*(x*x + y*y))); //elevation
    float t0 = 2.0*(w*y - z*x);
    if(t0 > 1.0)
        t0 = 1;
    else if(t0 < -1.0)
        t0 = -1;

    y = asin(t0);
    z = atan2((2.0*(w*z + y*y)), (1.0 - 2.0*(y*y + z*z))); //azimuth

    vector<float> ang {x, y, z};
    return ang;
}


vector<float> process_epoch(MatrixXd epoch_data, Compactor compactor, MatrixXd odo_csv, ofstream &logfile)
{
    float pi2 = 2*M_PI;
    Epoch epoch(compactor.dataset);
    epoch.dataset = compactor.dataset;
    epoch.data = epoch_data;
    epoch.ts = compactor.curr_ts;
    epoch.get_odo_data(odo_csv, compactor.dataset);

    if(compactor.radar_orients.size() == 0)
        compactor.radar_orients = get_radar_orient();

    //get velocity of current epoch relative to the moving frame
    if(compactor.test_suite.run_ls)
    {
        //run the naive ls computation
        vector<float> head_velo_ls = get_epoch_velo(compactor, epoch.data, true);
        float epoch_heading_ls = head_velo_ls[0];
        float x_velo_ls =        head_velo_ls[1];
        float y_velo_ls =        head_velo_ls[2];

        if(epoch_heading_ls < 0)
            epoch_heading_ls += 2*M_PI;

        float a = compactor.ls_heading + epoch_heading_ls;
        compactor.ls_heading = a - floor(a/pi2)*pi2;

        float speed_ls = pow(abs(x_velo_ls*x_velo_ls + y_velo_ls*y_velo_ls), 0.5);
        compactor.new_pos(speed_ls, true);

        compactor.heading = compactor.ls_heading;
        compactor.x_pos = compactor.ls_x_pos;
        compactor.y_pos = compactor.ls_y_pos;

        epoch.heading = epoch_heading_ls;
        epoch.x_velo = x_velo_ls;
        epoch.y_velo = y_velo_ls;
    }

    else
    {
        vector<float> head_velo = get_epoch_velo(compactor, epoch.data);
        epoch.heading = head_velo[0];
        epoch.x_velo =  head_velo[1];
        epoch.y_velo =  head_velo[2];

        if(epoch.heading < 0)
            epoch.heading += 2*M_PI;

        float a = compactor.heading + epoch.heading;
        compactor.ls_heading = a - floor(a/pi2)*pi2;

        float speed = pow(abs(epoch.x_velo*epoch.x_velo + epoch.y_velo*epoch.y_velo), 0.5);
        compactor.new_pos(speed);
    }

    compactor.log_epoch(epoch, logfile);
    vector<float> nh = new_pos_nh(compactor, epoch.x_velo, epoch.y_velo);
    return nh;
}


vector<float> new_pos_nh(Compactor compactor, float v_x, float v_y)
{
    float delta_t = compactor.next_ts - compactor.curr_ts;
    float y_pos_nh = -v_y*delta_t;

    float x_pos_nh;
    if(compactor.dataset == 2)
        x_pos_nh = v_x*delta_t;
    else
        x_pos_nh = -v_x*delta_t;

    vector<float> out = {x_pos_nh, y_pos_nh};
    return out;
}


vector<float> get_radar_orient(vector<float> cal_adjusts, int radar_idx, float adjust)
{
    vector<float> radar_orient = {0,0,0,0,0,0};
    radar_orient[0] = cal_adjusts[0] + 0;
    radar_orient[1] = cal_adjusts[1] + 0;
    radar_orient[2] = cal_adjusts[2] + M_PI;
    radar_orient[3] = cal_adjusts[3] + M_PI;
    radar_orient[4] = cal_adjusts[4] + 1.0471975493043784;
    radar_orient[5] = cal_adjusts[5] - 0.8726652191830278;

    if(radar_idx != -1)
    {
        radar_orient[radar_idx] += adjust;
    }

    return radar_orient;
}


vector<float> get_epoch_velo(Compactor compactor, MatrixXd epoch_data, bool ls)
{
    vector<vector<vector<double>>> radar_dict = get_radar_dict(epoch_data);

    MatrixXd A_SR = create_matrices(compactor, radar_dict);
    MatrixXd A  = A_SR.leftCols(2);
    VectorXd speed_radials = A_SR.col(2);

    if(ls)
        return least_squares_naive(A, speed_radials);
    else
        return least_squares_ransac(A, speed_radials);
}


vector<vector<vector<double>>> get_radar_dict(MatrixXd epoch_data)
{
    vector<vector<vector<double>>> radar_dict {{},{},{},{},{},{}};

    for(int r=0; r<epoch_data.rows(); r++)
    {
        VectorXd line = epoch_data.row(r);
        int radar_idx = line[9];
        radar_dict[radar_idx].push_back({});

        for(int c=0; c<line.size(); c++)
            radar_dict[radar_idx][ radar_dict[radar_idx].size()-1 ].push_back(line[c]);
    }
    return radar_dict;
}


MatrixXd create_matrices(Compactor compactor, vector<vector<vector<double>>> radar_dict)
{
    MatrixXd A_SR (0,3);

    for(int i=0; i<6; i++)
    {
        vector<vector<double>> curr_radar_data = radar_dict[i];
        float curr_quat_orient = compactor.radar_orients[i];

        for(int j=0; j<curr_radar_data.size(); j++)
        {
            vector<double> line = curr_radar_data[j];
            double theta_og = line[11] + curr_quat_orient;

            Vector3d new_row;
            new_row << cos(theta_og), sin(theta_og), line[6];

            A_SR.conservativeResize(A_SR.rows()+1, A_SR.cols());
            A_SR.row(A_SR.rows()-1) = new_row;
        }
    }
    return A_SR;
}


vector<float> least_squares_naive(MatrixXd A, VectorXd speed_radials)
{
    MatrixXd least_squares_velo = (A.transpose()*A).inverse() * A.transpose()*speed_radials;
    float ls_heading = atan2(least_squares_velo(1), least_squares_velo(0));
    vector<float> head_velo {ls_heading, least_squares_velo(0), least_squares_velo(1)};

    return head_velo;
}


vector<float> least_squares_ransac(MatrixXd A, VectorXd speed_radials)
{
    return least_squares_naive(A, speed_radials);
}
