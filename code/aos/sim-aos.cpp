#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
using namespace std;

struct Body{
    /*class attributes*/
    double position_x;
    double position_y;
    double position_z;
    double velocity_x;
    double velocity_y;
    double velocity_z;
    double mass;
    /*class methods*/
    Body(double position_x, double position_y, double position_z, double velocity_x, double velocity_y,
         double velocity_z, double mass) {
        //constructor (random for each body)
        Body::position_x = position_x;
        Body::position_y = position_y;
        Body::position_z = position_z;
        Body::velocity_x = velocity_x;
        Body::velocity_y = velocity_y;
        Body::velocity_z = velocity_z;
        Body::mass = mass;
    }
    ~Body() = default;
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
    mt19937_64 gen(random_seed);
    uniform_real_distribution<double> urd(min_pos_value, size_enclosure);
    normal_distribution<double> nd(1e+21,1e+15);
    vector<Body> vector_bodies;
    vector_bodies.reserve(num_objects);
    for (auto i = 0; i < num_objects; ++i){
        double pos_x = urd(gen);
        double pos_y = urd(gen);
        double pos_z = urd(gen);
        double mass = nd(gen);
        Body body(pos_x, pos_y, pos_z, 0.0, 0.0, 0.0, mass);
        vector_bodies.emplace_back(body);
    }
    //WRITING INITIAL CONFIGURATION-------------------------------------------------------------------------------------
    string init_config = "init_config_aos.txt";
    ofstream init_file(init_config);
    init_file << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << vector_bodies.size() << endl;
    for(auto &i: vector_bodies){
        init_file << fixed << setprecision(3) << i.position_x << " " << i.position_y << " " << i.position_z
                  << " " << i.velocity_x << " " << i.velocity_y << " " << i.velocity_z << " " << i.mass << endl;
    }
    //CHECKING FOR COLLISIONS INITIALLY---------------------------------------------------------------------------------
    for (unsigned long i = 0; i < vector_bodies.size()-1; ++i) {
        for (unsigned long j = i + 1; j < vector_bodies.size(); ++j) {
            if (vector_bodies[j].mass > 0) {
                double distance = std::sqrt((vector_bodies[i].position_x - vector_bodies[j].position_x) *
                                            (vector_bodies[i].position_x - vector_bodies[j].position_x) +
                                            (vector_bodies[i].position_y - vector_bodies[j].position_y) *
                                            (vector_bodies[i].position_y - vector_bodies[j].position_y) +
                                            (vector_bodies[i].position_z - vector_bodies[j].position_z) *
                                            (vector_bodies[i].position_z - vector_bodies[j].position_z));
                if (distance < 1) {
                    vector_bodies[i].mass += vector_bodies[j].mass;
                    vector_bodies[i].velocity_x += vector_bodies[j].velocity_x;
                    vector_bodies[i].velocity_y += vector_bodies[j].velocity_y;
                    vector_bodies[i].velocity_z += vector_bodies[j].velocity_z;
                    vector_bodies[j].mass = -1.0;
                }
            }
        }
    }
    vector_bodies.erase(std::remove_if(vector_bodies.begin(), vector_bodies.end(),
                                       [](Body &i) { return (i.mass < 0); }), vector_bodies.end());
    vector_bodies.shrink_to_fit();
    //ITERATIONS--------------------------------------------------------------------------------------------------------
    for (int iter = 0; iter < num_iterations; ++iter) {
        //CALCULATIONS OF FORCES, VELOCITIES, POSITIONS & REBOUNDS------------------------------------------------------
        vector<vector<double>> force(vector_bodies.size(), vector<double>(3, 0.0));
        for (unsigned long i = 0; i < vector_bodies.size(); ++i) {
            for (unsigned long j = i + 1; j < vector_bodies.size(); ++j) {
                if (i < vector_bodies.size() - 1){
                    double norm = std::sqrt((vector_bodies[j].position_x - vector_bodies[i].position_x) *
                                            (vector_bodies[j].position_x - vector_bodies[i].position_x) +
                                            (vector_bodies[j].position_y - vector_bodies[i].position_y) *
                                            (vector_bodies[j].position_y - vector_bodies[i].position_y) +
                                            (vector_bodies[j].position_z - vector_bodies[i].position_z) *
                                            (vector_bodies[j].position_z - vector_bodies[i].position_z));
                    force[i][0] += G * vector_bodies[i].mass * vector_bodies[j].mass * (
                            vector_bodies[j].position_x - vector_bodies[i].position_x) / norm / norm / norm;
                    force[i][1] += G * vector_bodies[i].mass * vector_bodies[j].mass * (
                            vector_bodies[j].position_y - vector_bodies[i].position_y) / norm / norm / norm;
                    force[i][2] += G * vector_bodies[i].mass * vector_bodies[j].mass * (
                            vector_bodies[j].position_z - vector_bodies[i].position_z) / norm / norm / norm;
                    /* fji = -fij */
                    force[j][0] -= G * vector_bodies[i].mass * vector_bodies[j].mass * (
                            vector_bodies[j].position_x - vector_bodies[i].position_x) / norm / norm / norm;
                    force[j][1] -= G * vector_bodies[i].mass * vector_bodies[j].mass * (
                            vector_bodies[j].position_y - vector_bodies[i].position_y) / norm / norm / norm;
                    force[j][2] -= G * vector_bodies[i].mass * vector_bodies[j].mass * (
                            vector_bodies[j].position_z - vector_bodies[i].position_z) / norm / norm / norm;
                }
            }
            /* update positions & velocities */
            vector_bodies[i].velocity_x += 1 / vector_bodies[i].mass * force[i][0] * time_step;
            vector_bodies[i].velocity_y += 1 / vector_bodies[i].mass * force[i][1] * time_step;
            vector_bodies[i].velocity_z += 1 / vector_bodies[i].mass * force[i][2] * time_step;
            vector_bodies[i].position_x += vector_bodies[i].velocity_x * time_step;
            vector_bodies[i].position_y += vector_bodies[i].velocity_y * time_step;
            vector_bodies[i].position_z += vector_bodies[i].velocity_z * time_step;
            //REBOUNDS--------------------------------------------------------------------------------------------------
            if (vector_bodies[i].position_x >= size_enclosure){
                vector_bodies[i].position_x = size_enclosure;
                vector_bodies[i].velocity_x = -vector_bodies[i].velocity_x;
            }
            else if (vector_bodies[i].position_x <= min_pos_value){
                vector_bodies[i].position_x = min_pos_value;
                vector_bodies[i].velocity_x = -vector_bodies[i].velocity_x;
            }
            if (vector_bodies[i].position_y >= size_enclosure){
                vector_bodies[i].position_y = size_enclosure;
                vector_bodies[i].velocity_y = -vector_bodies[i].velocity_y;
            }
            else if (vector_bodies[i].position_y <= min_pos_value){
                vector_bodies[i].position_y = min_pos_value;
                vector_bodies[i].velocity_y = -vector_bodies[i].velocity_y;
            }
            if (vector_bodies[i].position_z >= size_enclosure){
                vector_bodies[i].position_z = size_enclosure;
                vector_bodies[i].velocity_z = -vector_bodies[i].velocity_z;
            }
            else if (vector_bodies[i].position_z <= min_pos_value){
                vector_bodies[i].position_z = min_pos_value;
                vector_bodies[i].velocity_z = -vector_bodies[i].velocity_z;
            }
        }
        //COLLISIONS
        for (unsigned long i = 0; i < vector_bodies.size() - 1; ++i) {
            for (unsigned long j = i + 1; j < vector_bodies.size(); ++j) {
                if (vector_bodies[j].mass > 0) {
                    double distance = std::sqrt((vector_bodies[i].position_x - vector_bodies[j].position_x) *
                                                (vector_bodies[i].position_x - vector_bodies[j].position_x) +
                                                (vector_bodies[i].position_y - vector_bodies[j].position_y) *
                                                (vector_bodies[i].position_y - vector_bodies[j].position_y) +
                                                (vector_bodies[i].position_z - vector_bodies[j].position_z) *
                                                (vector_bodies[i].position_z - vector_bodies[j].position_z));
                    if (distance < 1) {
                        vector_bodies[i].mass += vector_bodies[j].mass;
                        vector_bodies[i].velocity_x += vector_bodies[j].velocity_x;
                        vector_bodies[i].velocity_y += vector_bodies[j].velocity_y;
                        vector_bodies[i].velocity_z += vector_bodies[j].velocity_z;
                        vector_bodies[j].mass = -1.0;
                    }
                }
            }
        }
        vector_bodies.erase(std::remove_if(vector_bodies.begin(), vector_bodies.end(),
                                           [](Body &i) { return (i.mass < 0); }), vector_bodies.end());
        vector_bodies.shrink_to_fit();
    }
    //WRITING FINAL CONFIGURATION---------------------------------------------------------------------------------------
    string final_config = "final_config_aos.txt";
    ofstream final_file(final_config);
    final_file << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << vector_bodies.size() << endl;
    for(auto &i: vector_bodies){
        final_file << fixed << setprecision(3) << i.position_x << " " << i.position_y << " " << i.position_z
                   << " " << i.velocity_x << " " << i.velocity_y << " " << i.velocity_z << " " << i.mass << endl;
    }
    return 0;
}