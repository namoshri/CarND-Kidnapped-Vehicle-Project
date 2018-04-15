/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::DumpParticle(struct Particle particle) {
	cout << "Particle " << particle.id + 1 << " " << particle.x << " " << particle.y <<  " " << particle.theta << " " << particle.weight << endl;
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.

	//"gen" is the random engine.
	default_random_engine gen;

	//creates a normal (Gaussian) distribution
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	//set num of particles and update each particle
	num_particles = 30;
	particles.resize(num_particles);
	for (int i = 0; i < num_particles; ++i) {
		particles[i].id = i;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
		particles[i].weight = 1;

		// Print your particle to the terminal.
		//DumpParticle(particles[i]);
	}
	//Inform that intialization done
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// adding noise through std::normal_distribution and std::default_random_engine useful.
	// Reference:
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//"gen" is the random engine.
	default_random_engine gen;

	//Temp variables for code simplicity and avoid multiple calc
	double yaw_rate_delta = yaw_rate * delta_t;
	double vel_yaw_rate = velocity / yaw_rate;

	for (int i = 0; i < num_particles; ++i) {
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		double theta_yaw_rate = theta + yaw_rate_delta;


		//updating x, y and the yaw angle when the yaw rate is not equal to zero:
		if (fabs(yaw_rate) > 0.0001) {
			//First get x and y
			x = x + (vel_yaw_rate * (sin(theta_yaw_rate) - sin(theta)));
			y = y + (vel_yaw_rate * (cos(theta) - cos(theta_yaw_rate)));

			//..and then update theta as theta is needed in x and y
			theta = theta + yaw_rate_delta;
		} else { //updating x, y as yaw rate is zero:
			x = x + (velocity * delta_t * cos(theta));
			y = y + (velocity * delta_t * sin(theta));
		}

		// normal distribution and update for predecition
		normal_distribution<double> dist_x(x, std_pos[0]);
		normal_distribution<double> dist_y(y, std_pos[1]);
		normal_distribution<double> dist_theta(theta, std_pos[2]);

		// Particle before prediction.
		//DumpParticle(particles[i]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

		// Particle ater prediction.
		//DumpParticle(particles[i]);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
