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
	weights.resize(num_particles, 1.0);

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
		if (fabs(yaw_rate) > 0.00001) {
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

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> preds, std::vector<LandmarkObs>& obs) {
	// Find the predicted measurement that is closest to each observed measurement and assign the 
	// observed measurement to this particular landmark.
	// used it as a helper during the updateWeights phase.

	for (unsigned int i = 0; i < obs.size(); i++) {
		double mind = 1.0e20; //Some high number
		for (unsigned int j = 0; j < preds.size(); j++) {
			double d = dist(obs[i].x, obs[i].y, preds[j].x, preds[j].y);
			if (d < mind) {
				obs[i].id = preds[j].id;
				mind = d;
			}
		}
	}
}

void ParticleFilter::predictList(struct Particle p, double sensor_range,
				const Map &maps,
				std::vector<LandmarkObs> *pred_landm) {

	std::vector<LandmarkObs> t_predicts;
	for (unsigned int i = 0; i < maps.landmark_list.size(); i++) {
		auto m = maps.landmark_list[i];
		double d = dist(p.x, p.y, m.x_f, m.y_f);
		
		if (d < sensor_range) {
			LandmarkObs temp;

			temp.id = m.id_i;
			temp.x = m.x_f;
			temp.y = m.y_f;
			t_predicts.push_back(temp);
		}
	}
	*pred_landm = t_predicts;
}

void ParticleFilter::transform(struct Particle p, const std::vector<LandmarkObs> &obs,
			std::vector<LandmarkObs> *trans_obs) {

//x_map= x_part + (np.cos(theta) * x_obs) - (np.sin(theta) * y_obs)
//y_map= y_part + (np.sin(theta) * x_obs) + (np.cos(theta) * y_obs)

	vector<LandmarkObs> t_transform;
	for (unsigned int i = 0; i < obs.size(); i++) {
		LandmarkObs temp;
		temp.x = p.x + (cos(p.theta) * obs[i].x) - (sin(p.theta) * obs[i].y);
		temp.y = p.y + (sin(p.theta) * obs[i].x) + (cos(p.theta) * obs[i].y);
		temp.id = obs[i].id;
		t_transform.push_back(temp);
	}
	*trans_obs = t_transform;
}

double ParticleFilter::particleWeight(std::vector<LandmarkObs> preds, std::vector<LandmarkObs>& obs,
			double constpart, double std_x, double std_y) {
	double finalp = 1.0;
	double tempp = 1.0;
	double stdxsq = 2.0 * std_x * std_x;
	double stdysq = 2.0 * std_y * std_y;

	for (unsigned int j = 0; j < obs.size(); j++) {
		LandmarkObs obsw = obs[j];
		LandmarkObs temp;
		double x = obsw.x;
		double y = obsw.y;

		for (unsigned int k = 0; k < preds.size(); k++) {
			if (obsw.id == preds[k].id) {
				temp = preds[k];
				break;
			}
		}
		auto dx = (temp.x - x) * (temp.x - x);
		auto dy = (temp.y - y) * (temp.y - y);
		tempp = exp (-(dx/stdxsq + dy/stdysq)) * constpart;
		finalp *= tempp;
	}
	return finalp;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	//   Update the weights of each particle using a mult-variate Gaussian distribution.
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	//   1. NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	//   Reference from past students on discussions.udacity.com to get flow and hints

	weights.clear();

	double std_x = std_landmark[0];
	double std_y = std_landmark[1];
	double constpart = (1.0 / (2.0 * M_PI * std_x * std_y));

	for (int i = 0; i < num_particles; i++) {
		//1. Transform
		std::vector<LandmarkObs> trans_obs;
		transform(particles[i], observations, &trans_obs);

		//2. list possible landmarks
		std::vector<LandmarkObs> pred_lm;
		predictList(particles[i], sensor_range, map_landmarks, &pred_lm);

		//3. data association
		dataAssociation(pred_lm, trans_obs);

		//4. Weight calcualtion
		particles[i].weight = particleWeight(pred_lm, trans_obs, constpart, std_x, std_y);
		weights.push_back(particles[i].weight);
	}
}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight. 
	// References : std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// Refered past student from discussion.udacity

	//"gen" is the random engine.
	default_random_engine gen;
	std::discrete_distribution<> dd(weights.begin(), weights.end());
	std::vector<Particle> temp(num_particles);

	for (int i = 0; i < num_particles; i++)
		temp[i] = particles[dd(gen)];
	particles = temp;
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
