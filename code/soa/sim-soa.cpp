#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
using namespace std;

struct Bodies {
    /*class attributes*/
    vector<double> positions_x;
    vector<double> positions_y;
    vector<double> positions_z;
    vector<double> velocities_x;
    vector<double> velocities_y;
    vector<double> velocities_z;
    vector<double> masses;
    /*class methods*/
    explicit Bodies(int num_bodies) : positions_x(num_bodies), positions_y(num_bodies), positions_z(num_bodies),
    velocities_x(num_bodies, 0.0), velocities_y(num_bodies, 0.0), velocities_z(num_bodies, 0.0),
    masses(num_bodies){};
    ~Bodies() = default;
};

int main(int argc, char* argv[]){
    // CHECKING ERRORS -------------------------------------------------------------------------------------------------
    // first, we check whether the number of arguments is indeed 5
    string arg_names[5] = {"num_objects", "num_iterations", "random_seed", "size_enclosure", "time_step"};
    if (argc != 6){
        cerr << "The number of arguments is incorrect. You have introduced the next "
             << (argc -1) << " arguments:" << endl;
        if (argv[1] != nullptr){
            cerr << arg_names[0] << ": " << argv[1] << endl;
        } else {
            cerr << arg_names[0] << ": " << "??" << endl;
        }
        if (argv[1] != nullptr && argv[2] != nullptr){
            cerr << arg_names[1] << ": " << argv[2] << endl;
        } else {
            cerr << arg_names[1] << ": " << "??" << endl;
        }
        if (argv[1] != nullptr && argv[2] != nullptr && argv[3] != nullptr){
            cerr << arg_names[2] << ": " << argv[3] << endl;
        } else {
            cerr << arg_names[2] << ": " << "??" << endl;
        }
        if (argv[1] != nullptr && argv[2] != nullptr && argv[3] != nullptr
            && argv[4] != nullptr){
            cerr << arg_names[3] << ": " << argv[4] << endl;
        } else {
            cerr << arg_names[3] << ": " << "??" << endl;
        }
        if (argv[1] != nullptr && argv[2] != nullptr && argv[3] != nullptr
            && argv[4] != nullptr && argv[5] != nullptr){
            cerr << arg_names[4] << ": " << argv[5] << endl;
        } else {
            cerr << arg_names[4] << ": " << "??" << endl;
        }
        cerr << "You have to introduce 5 arguments. Please, introduce the missing ones\n";
        return -1;
    }
    /* now we check whether the format of the arguments is correct;
     * to do so, we first define the default error handling message that we will be using*/
    const char *error_handler = "\nWrong input parameter(s). The input parameters must have the following format and order:\n"
                                "num_objects: positive integer number greater than 0\n"
                                "num_iterations: positive integer number greater than 0\n"
                                "random_seed: positive integer number greater than 0\n"
                                "size_enclosure: positive real number greater than 0\n"
                                "time_step: positive real number greater than 0\n";
    // here we define the wrong arguments message
    const char *wrong_arg_half_1 = "Please, introduce again the ";
    const char *wrong_arg_half_2 = "parameter(s)";

    // first, we save all the arguments into variables
    vector<char *> p(5, nullptr);
    vector<string> s = {argv[1], argv[2], argv[3], argv[4], argv[5]};
    vector<string> wrong_args;
    const int num_objects = (int)strtol(s[0].data(), &p[0], 10);
    const int num_iterations = (int)strtol(s[1].data(), &p[1], 10);
    const int random_seed = (int)strtol(s[2].data(), &p[2], 10);
    const auto size_enclosure = (double)strtof(s[3].data(), &p[3]);
    const auto time_step = (double)strtof(s[4].data(), &p[4]);

    // now, we can check whether each of them is correct
    if (num_objects <= 0 || (typeid(num_objects) != typeid(int))|| p[0] != (s[0].data() + s[0].size())){
        wrong_args.emplace_back(arg_names[0]);
    }
    if (num_iterations < 0 || (typeid(num_iterations) != typeid(int)) || p[1] != (s[1].data() + s[1].size())){
        wrong_args.emplace_back(arg_names[1]);
    }
    if (random_seed < 0 || (typeid(random_seed) != typeid(int)) || p[2] != (s[2].data() + s[2].size())){
        wrong_args.emplace_back(arg_names[2]);
    }
    if (size_enclosure < 0 || (typeid(size_enclosure) != typeid(double)) || p[3] != (s[3].data() + s[3].size())){
        wrong_args.emplace_back(arg_names[3]);
    }
    if (time_step< 0 || (typeid(time_step) != typeid(double)) || p[4] != (s[4].data() + s[4].size())){
        wrong_args.emplace_back(arg_names[4]);
    }
    // print wrong arguments message
    if (!wrong_args.empty()){
        cerr << error_handler << endl;
        cerr << wrong_arg_half_1;
        for (auto & wrong_arg : wrong_args) {
            cerr << '"' << wrong_arg << '"' << " ";
        }
        cerr << wrong_arg_half_2 << endl;
        return -2;
    }
    //PRINTING INPUT ARGUMENTS------------------------------------------------------------------------------------------
    cout << "Creating simulation:" << "\n  " << "num_objects: " << num_objects << "\n  "
         << "num_iterations: " << num_iterations << "\n  " << "random_seed: " << random_seed << "\n  "
         << "size_enclosure: " << size_enclosure << "\n  " << "time_step: " << time_step << endl;
    //CREATING SIMULATION PARAMETERS------------------------------------------------------------------------------------
    const double G = 6.674e-11;
    const double min_pos_value = 0.0;
    Bodies bodies(num_objects);
    mt19937_64 gen(random_seed);
    uniform_real_distribution<double> urd(min_pos_value, size_enclosure);
    normal_distribution<double> nd(1e+21,1e+15);
    for (auto i = 0; i < num_objects; ++i){
        bodies.positions_x[i] = urd(gen);
        bodies.positions_y[i] = urd(gen);
        bodies.positions_z[i] = urd(gen);
        bodies.masses[i] = nd(gen);
    }
    //WRITING INITIAL CONFIGURATION-------------------------------------------------------------------------------------
    string init_config = "init_config_soa.txt";
    ofstream init_file(init_config);
    init_file << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << bodies.masses.size() << endl;
    for(unsigned long i = 0; i < bodies.masses.size(); ++i){
        init_file << fixed << setprecision(3) << bodies.positions_x[i] << " " << bodies.positions_y[i] << " "
        << bodies.positions_z[i] << " " << bodies.velocities_x[i] << " " << bodies.velocities_y[i] << " "
        << bodies.velocities_z[i] << " " << bodies.masses[i] << endl;
    }
    //CHECKING FOR COLLISIONS INITIALLY---------------------------------------------------------------------------------
    for (unsigned long i = 0; i < bodies.masses.size()-1; ++i) {
        for (unsigned long j = i + 1; j < bodies.masses.size(); ++j) {
            if (bodies.masses[j] > 0) {
                double distance = std::sqrt((bodies.positions_x[i] - bodies.positions_x[j]) *
                                            (bodies.positions_x[i] - bodies.positions_x[j]) +
                                            (bodies.positions_y[i] - bodies.positions_y[j]) *
                                            (bodies.positions_y[i] - bodies.positions_y[j]) +
                                            (bodies.positions_z[i] - bodies.positions_z[j]) *
                                            (bodies.positions_z[i] - bodies.positions_z[j]));
                if (distance < 1) {
                    bodies.masses[i] += bodies.masses[j];
                    bodies.velocities_x[i] += bodies.velocities_x[j];
                    bodies.velocities_y[i] += bodies.velocities_y[j];
                    bodies.velocities_z[i] += bodies.velocities_z[j];
                    bodies.masses[j] = -1.0;
                    bodies.positions_x[j] = -1.0;
                    bodies.positions_y[j] = -1.0;
                    bodies.positions_z[j] = -1.0;
                    bodies.velocities_x[j] = -1.0;
                    bodies.velocities_y[j] = -1.0;
                    bodies.velocities_z[j] = -1.0;
                }
            }
        }
    }
    bodies.masses.erase(std::remove_if(bodies.masses.begin(), bodies.masses.end(),
                                       [](double &i) { return (i < 0); }), bodies.masses.end());
    bodies.masses.shrink_to_fit();

    bodies.positions_x.erase(std::remove_if(bodies.positions_x.begin(), bodies.positions_x.end(),
                                       [](double &i) { return (i < 0); }), bodies.positions_x.end());
    bodies.positions_x.shrink_to_fit();
    bodies.positions_y.erase(std::remove_if(bodies.positions_y.begin(), bodies.positions_y.end(),
                                       [](double &i) { return (i < 0); }), bodies.positions_y.end());
    bodies.positions_y.shrink_to_fit();
    bodies.positions_z.erase(std::remove_if(bodies.positions_z.begin(), bodies.positions_z.end(),
                                       [](double &i) { return (i < 0); }), bodies.positions_z.end());
    bodies.positions_z.shrink_to_fit();

    bodies.velocities_x.erase(std::remove_if(bodies.velocities_x.begin(), bodies.velocities_x.end(),
                                            [](double &i) { return (i < 0); }), bodies.velocities_x.end());
    bodies.velocities_x.shrink_to_fit();
    bodies.velocities_y.erase(std::remove_if(bodies.velocities_y.begin(), bodies.velocities_y.end(),
                                            [](double &i) { return (i < 0); }), bodies.velocities_y.end());
    bodies.velocities_y.shrink_to_fit();
    bodies.velocities_z.erase(std::remove_if(bodies.velocities_z.begin(), bodies.velocities_z.end(),
                                            [](double &i) { return (i < 0); }), bodies.velocities_z.end());
    bodies.velocities_z.shrink_to_fit();
    //ITERATIONS--------------------------------------------------------------------------------------------------------
    for (int iter = 0; iter < num_iterations; ++iter) {
        //CALCULATIONS OF FORCES, VELOCITIES, POSITIONS & REBOUNDS------------------------------------------------------
        vector<vector<double>> force(bodies.masses.size(), vector<double>(3, 0.0));
        for (unsigned long i = 0; i < bodies.masses.size(); ++i) {
            for (unsigned long j = i + 1; j < bodies.masses.size(); ++j) {
                if (i < bodies.masses.size() - 1){
                    double norm = std::sqrt((bodies.positions_x[j] - bodies.positions_x[i]) *
                                            (bodies.positions_x[j] - bodies.positions_x[i]) +
                                            (bodies.positions_y[j] - bodies.positions_y[i]) *
                                            (bodies.positions_y[j] - bodies.positions_y[i]) +
                                            (bodies.positions_z[j] - bodies.positions_z[i]) *
                                            (bodies.positions_z[j] - bodies.positions_z[i]));
                    force[i][0] += G * bodies.masses[i] * bodies.masses[j] * (bodies.positions_x[j] - bodies.positions_x[i]) / norm / norm / norm;
                    force[i][1] += G * bodies.masses[i] * bodies.masses[j] * (bodies.positions_y[j] - bodies.positions_y[i]) / norm / norm / norm;
                    force[i][2] += G * bodies.masses[i] * bodies.masses[j] * (bodies.positions_z[j] - bodies.positions_z[i]) / norm / norm / norm;
                    /* fji = -fij */
                    force[j][0] -= G * bodies.masses[i] * bodies.masses[j] * (bodies.positions_x[j] - bodies.positions_x[i]) / norm / norm / norm;
                    force[j][1] -= G * bodies.masses[i] * bodies.masses[j] * (bodies.positions_y[j] - bodies.positions_y[i]) / norm / norm / norm;
                    force[j][2] -= G * bodies.masses[i] * bodies.masses[j] * (bodies.positions_z[j] - bodies.positions_z[i]) / norm / norm / norm;
                }
            }
            /* update positions & velocities */
            bodies.velocities_x[i] += 1 / bodies.masses[i] * force[i][0] * time_step;
            bodies.velocities_y[i] += 1 / bodies.masses[i] * force[i][1] * time_step;
            bodies.velocities_z[i] += 1 / bodies.masses[i] * force[i][2] * time_step;
            bodies.positions_x[i] += bodies.velocities_x[i] * time_step;
            bodies.positions_y[i] += bodies.velocities_y[i] * time_step;
            bodies.positions_z[i] += bodies.velocities_z[i] * time_step;
            //REBOUNDS--------------------------------------------------------------------------------------------------
            if (bodies.positions_x[i] >= size_enclosure){
                bodies.positions_x[i] = size_enclosure;
                bodies.velocities_x[i] = -bodies.velocities_x[i];
            }
            else if (bodies.positions_x[i] <= min_pos_value){
                bodies.positions_x[i] = min_pos_value;
                bodies.velocities_x[i] = -bodies.velocities_x[i];
            }
            if (bodies.positions_y[i] >= size_enclosure){
                bodies.positions_y[i] = size_enclosure;
                bodies.velocities_y[i] = -bodies.velocities_y[i];
            }
            else if (bodies.positions_y[i] <= min_pos_value){
                bodies.positions_y[i] = min_pos_value;
                bodies.velocities_y[i] = -bodies.velocities_y[i];
            }
            if (bodies.positions_z[i] >= size_enclosure){
                bodies.positions_z[i] = size_enclosure;
                bodies.velocities_z[i] = -bodies.velocities_z[i];
            }
            else if (bodies.positions_z[i] <= min_pos_value){
                bodies.positions_z[i] = min_pos_value;
                bodies.velocities_z[i] = -bodies.velocities_z[i];
            }
        }
        //COLLISIONS
        for (unsigned long i = 0; i < bodies.masses.size()-1; ++i) {
            for (unsigned long j = i + 1; j < bodies.masses.size(); ++j) {
                if (bodies.masses[j] > 0) {
                    double distance = std::sqrt((bodies.positions_x[i] - bodies.positions_x[j]) *
                                                (bodies.positions_x[i] - bodies.positions_x[j]) +
                                                (bodies.positions_y[i] - bodies.positions_y[j]) *
                                                (bodies.positions_y[i] - bodies.positions_y[j]) +
                                                (bodies.positions_z[i] - bodies.positions_z[j]) *
                                                (bodies.positions_z[i] - bodies.positions_z[j]));
                    if (distance < 1) {
                        bodies.masses[i] += bodies.masses[j];
                        bodies.velocities_x[i] += bodies.velocities_x[j];
                        bodies.velocities_y[i] += bodies.velocities_y[j];
                        bodies.velocities_z[i] += bodies.velocities_z[j];
                        bodies.masses[j] = -1.0;
                        bodies.positions_x[j] = -1.0;
                        bodies.positions_y[j] = -1.0;
                        bodies.positions_z[j] = -1.0;
                        bodies.velocities_x[j] = 0.0;
                        bodies.velocities_y[j] = 0.0;
                        bodies.velocities_z[j] = 0.0;
                    }
                }
            }
        }
        bodies.masses.erase(std::remove_if(bodies.masses.begin(), bodies.masses.end(),
                                           [](double &i) { return (i < 0); }), bodies.masses.end());
        bodies.masses.shrink_to_fit();

        bodies.positions_x.erase(std::remove_if(bodies.positions_x.begin(), bodies.positions_x.end(),
                                                [](double &i) { return (i < 0); }), bodies.positions_x.end());
        bodies.positions_x.shrink_to_fit();
        bodies.positions_y.erase(std::remove_if(bodies.positions_y.begin(), bodies.positions_y.end(),
                                                [](double &i) { return (i < 0); }), bodies.positions_y.end());
        bodies.positions_y.shrink_to_fit();
        bodies.positions_z.erase(std::remove_if(bodies.positions_z.begin(), bodies.positions_z.end(),
                                                [](double &i) { return (i < 0); }), bodies.positions_z.end());
        bodies.positions_z.shrink_to_fit();

        bodies.velocities_x.erase(std::remove_if(bodies.velocities_x.begin(), bodies.velocities_x.end(),
                                                 [](double &i) { return (i == 0); }), bodies.velocities_x.end());
        bodies.velocities_x.shrink_to_fit();
        bodies.velocities_y.erase(std::remove_if(bodies.velocities_y.begin(), bodies.velocities_y.end(),
                                                 [](double &i) { return (i == 0); }), bodies.velocities_y.end());
        bodies.velocities_y.shrink_to_fit();
        bodies.velocities_z.erase(std::remove_if(bodies.velocities_z.begin(), bodies.velocities_z.end(),
                                                 [](double &i) { return (i == 0); }), bodies.velocities_z.end());
        bodies.velocities_z.shrink_to_fit();
    }
    //WRITING FINAL CONFIGURATION---------------------------------------------------------------------------------------
    string final_config = "final_config_soa.txt";
    ofstream final_file(final_config);
    final_file << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << bodies.masses.size() << endl;
    for(unsigned long i = 0; i < bodies.masses.size(); ++i){
        final_file << fixed << setprecision(3) << bodies.positions_x[i] << " " << bodies.positions_y[i] << " "
                  << bodies.positions_z[i] << " " << bodies.velocities_x[i] << " " << bodies.velocities_y[i] << " "
                  << bodies.velocities_z[i] << " " << bodies.masses[i] << endl;
    }
    return 0;
}