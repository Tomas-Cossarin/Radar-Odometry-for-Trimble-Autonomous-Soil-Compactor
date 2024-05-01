#include "functions.h"
#include "classes.h"


Radar::Radar(double radar_pos_rax, double radar_pos_ray)
{
    x_pos_ra = radar_pos_rax;
    y_pos_ra = radar_pos_ray;
}

double Radar::get_theta_ra(double depth, double theta_og)
{
    double pos_ogx = depth*cos(theta_og);
    double pos_ogy = depth*sin(theta_og);
    double pos_rax = pos_ogx + x_pos_ra;
    double pos_ray = pos_ogy + y_pos_ra;
    double theta_ra = atan2(pos_ray, pos_rax);

    if(theta_ra < 0)
    {
        theta_ra += 2*M_PI;
    }
    if(pos_rax<0 && pos_ray>0)
    {
        theta_ra += M_PI/2;
    }
    else if(pos_rax<0 && pos_ray<0)
    {
        theta_ra += M_PI;
    }
    else if(pos_rax>0 && pos_ray<0)
    {
        theta_ra += 1.5*M_PI;
    }
    return theta_ra;
}


Epoch::Epoch(int dataset, float ts, float heading, float x_velo, float y_velo)
{
    this->ts = ts;
    this->heading = heading;
    this->x_velo = x_velo;
    this->y_velo = y_velo;
    this->dataset = dataset;
}

void Epoch::get_odo_data(MatrixXd odo_csv, int ds)
{
    // epoch must be set to valid timestep before this is called
    odo_data = get_odo_ts_data(ts, odo_csv);
    odo_x_pos = odo_data[1];
    odo_y_pos = odo_data[2];
    odo_heading = get_odo_heading_ts(odo_data, ds);
}


TestSuite::TestSuite()
{
    prep_ransac = false;
    ransac = false;
    calibration = false;
    no_ransac = false;
    run_ls = false;
}

void TestSuite::calibration_test(bool ransac)
{
    prep_ransac = false;
    ransac = false;
    calibration = true;
    no_ransac = true;
    run_ls = true;

    if (ransac == true)
    {
        no_ransac = false;
        run_ls = false;
    }
}

void TestSuite::ransac_test(bool create_outlier_data)
{
    prep_ransac = create_outlier_data;
    calibration = false;
    no_ransac = false;
    ransac = true;
    run_ls = true;
}

void TestSuite::ls_only()
{
    prep_ransac = false;
    ransac = false;
    calibration = false;
    no_ransac = true;
    run_ls = true;
}


Compactor::Compactor(ofstream &logfile, int dataset, float x_pos, float y_pos, float x_pos_nh, float y_pos_nh,
          float heading, float curr_ts, float next_ts, bool run_ls)
{
    // position is given in meters relative to the static local frame
    this->dataset = dataset;
    this->test_suite = TestSuite();

    this->x_pos = x_pos;
    this->y_pos = y_pos;
    this->heading = heading;
    this->x_pos_nh = x_pos_nh;
    this->y_pos_nh = y_pos_nh;
    this->epoch_start = 0;
    this->curr_ts = curr_ts;
    this->next_ts = next_ts;
    this->ls_heading = heading;
    this->ls_x_pos = x_pos;
    this->ls_y_pos = y_pos;
    this->run_ls = run_ls;
    this->get_radar_data_col_lst();

    string str = "timestep,duration,speed,x_velo_mf,y_velo_mf,x_pos_sf,y_pos_sf,x_pos_nh,y_pos_nh,heading_mf,heading_sf,x_pos_odo,y_pos_odo,heading_odo";
    write_to_log(logfile, str);
}

void Compactor::get_radar_data_col_lst()
{
    vector<string> col_lst;
    col_lst = {"source_name", "time", "id", "depth", "width", "height",
                   "speed_radial", "snr", "rcs", "radar_id", "og_range", "og_az", "og_el", "og_sr"};

    this->radar_data_col_lst = col_lst;
}

void Compactor::new_pos(float speed, bool ls)
{
    /*
    Calculates expected position based on current speed and heading

    params:
        old_coord:  old coordinates in (x, y) format
        heading:    direction facing in radians between -pi and pi
        speed:      speed in m/s
        t_step:     timestep in seconds

    Return:
        new_coord: expected coordinates in (x, y) format
    */
    float delta_t = this->next_ts - this->curr_ts;

    if (!ls)
    {
        assert(0 <= this->heading && this->heading <= 2*M_PI);
        this->x_pos += cos(this->heading) * speed * delta_t;
        this->y_pos += sin(this->heading) * speed * delta_t;
    }
    else
    {
        assert(0 <= this->ls_heading && this->ls_heading <= 2*M_PI);
        this->ls_x_pos += cos(this->ls_heading) * speed * delta_t;
        this->ls_y_pos += sin(this->ls_heading) * speed * delta_t;
    }
}

void Compactor::init_radar_positions()
{
    radar_lst.push_back(Radar({1.13, 0.97}));       //0: front left
    radar_lst.push_back(Radar({1.13, -0.97}));      //1: front right

    radar_lst.push_back(Radar({-1.981, 0.445}));    //2: rear left
    radar_lst.push_back(Radar({-1.981, -0.445}));   //3: rear right

    radar_lst.push_back(Radar({-0.902, 0.921}));    //4: side left
    radar_lst.push_back(Radar({-0.902, -0.921}));   //5: side right
}

void Compactor::log_epoch(Epoch epoch, ofstream &logfile)
{
    double speed = sqrt(abs(epoch.x_velo * epoch.x_velo + epoch.y_velo * epoch.y_velo));
    float delta_t = next_ts - curr_ts;

    write_to_log(logfile,
        to_string(curr_ts) + "," +
        to_string(delta_t) + "," +
        to_string(speed) + "," +
        to_string(epoch.x_velo) + "," +
        to_string(epoch.y_velo) + "," +
        to_string(x_pos) + "," +
        to_string(y_pos) + "," +
        to_string(x_pos_nh) + "," +
        to_string(y_pos_nh) + "," +
        to_string(epoch.heading) + "," +
        to_string(heading) + "," +
        to_string(epoch.odo_x_pos) + "," +
        to_string(epoch.odo_y_pos) + "," +
        to_string(epoch.odo_heading));
}
