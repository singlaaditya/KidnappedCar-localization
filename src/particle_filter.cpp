/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
  

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  std::default_random_engine gen;
  std::normal_distribution<double> dist_motion_x(0, std[0]);
  std::normal_distribution<double> dist_motion_y(0, std[1]);
  std::normal_distribution<double> dist_theta(0, std[2]);
  
  for (int i = 0; i < num_particles; i++) {
  	Particle p;
    double position_variance_x = dist_motion_x(gen);
    double position_variance_y = dist_motion_y(gen);
    double yaw_variance_comp = dist_theta(gen);
    p.id = i;
    p.x = x + position_variance_x;
    p.y = y + position_variance_y;
    p.theta = theta + yaw_variance_comp;
    p.weight = 1;
    particles.push_back(p);
  }
  
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  std::normal_distribution<double> dist_motion_x(0.0, std_pos[0]);
  std::normal_distribution<double> dist_motion_y(0.0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0.0, std_pos[2]);
  
  for (int i =0; i < particles.size(); i++) {
       if (fabs(yaw_rate) < 0.00001) {  
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    }
    else {
   particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
                                                                                                                      
   particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
                                              
   particles[i].theta += yaw_rate*delta_t;
    }		
    
    particles[i].x += dist_motion_x(gen);
    particles[i].y += dist_motion_y(gen);
    particles[i].theta += dist_theta(gen);
}}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    for (int i = 0; i < map_landmarks.landmark_list.size(); i++) {
    	std::cout << i << " " << map_landmarks.landmark_list[i].id_i << std::endl;
    }
    for (int i =0; i < particles.size(); i++) {
	    particles[i].weight = 1;
      
     vector<LandmarkObs> predictions;
    for(const auto& lm: map_landmarks.landmark_list){
      double distance = dist(particles[i].x, particles[i].y, lm.x_f, lm.y_f);
      if( distance < sensor_range){ // if the landmark is within the sensor range, save it to predictions
        predictions.push_back(LandmarkObs{lm.id_i, lm.x_f, lm.y_f});
      }
    }
      	vector<LandmarkObs> transformed;
      	std::vector<int> associations;
      	std::vector<double> sense_x;
      	std::vector<double> sense_y;
    	for (int j = 0; j < observations.size(); j++) {
        	LandmarkObs temp;
          	temp.id = observations[j].id;
          	temp.x = cos(particles[i].theta)*observations[j].x - sin(particles[i].theta) * observations[j].y + particles[i].x;
            temp.y = sin(particles[i].theta)*observations[j].x + cos(particles[i].theta) * observations[j].y + particles[i].y;
            transformed.push_back(temp);
        }
        for (int j = 0; j <transformed.size(); j++) {
        	double min_distance = std::numeric_limits<float>::max() ;
          	int map_id = -1;
            for (int k = 0; k < predictions.size(); k++){
            	double distance = dist(transformed[j].x, transformed[j].y, predictions[k].x, predictions[k].y);
             	if (distance < min_distance){
                 	min_distance = distance;
                  	map_id = predictions[k].id;
                }
            	
            }
            associations.push_back(map_id);
            sense_x.push_back(transformed[j].x);
          	sense_y.push_back(transformed[j].y);
          	
        }
      	SetAssociations(particles[i], associations, sense_x, sense_y);
            std::cout << "Aditya" << associations.size()  << sense_x.size() << sense_y.size() <<" " << predictions.size()<< std::endl;
          	for (int j = 0; j < associations.size(); j++) {
        	particles[i].weight *= exp(-1 * ( pow((sense_x[j] - map_landmarks.landmark_list[associations[j] - 1].x_f), 2)/(2 * std_landmark[0]* std_landmark[0]) + pow((sense_y[j] - map_landmarks.landmark_list[associations[j] - 1].y_f), 2)/(2 * std_landmark[1] * std_landmark[1]) )) / (2 * M_PI * std_landmark[0] * std_landmark[1]);  
        }
        
    }
  	double probability_sum = 0;
     for ( int i =0 ; i < particles.size(); i++) {
     	probability_sum += particles[i].weight;
     }
     for ( int i =0 ; i < particles.size(); i++) {
     	particles[i].weight /= probability_sum;
        weights.push_back(particles[i].weight);
     }
   	// Calculate weight based on associations
  /*
    for (int i =0; i < particles.size(); i++) {
        auto associations = particles[i].associations;
      	auto sense_x = particles[i].sense_x;
      	auto sense_y = particles[i].sense_y;
    	for (int j = 0; j < associations.size(); j++) {
        	particles[i].weight *= exp(-1 * ( pow((sense_x[j] - predictions[associations[j]].x), 2)/(2 * std_landmark[0]* std_landmark[0]) + pow((sense_y[j] - predictions[associations[j]].y), 2)/(2 * std_landmark[1] * std_landmark[1]) )) / (2 * M_PI * std_landmark[0] * std_landmark[1]);  
        }
        weights.push_back(particles[i].weight);
    } */
}

void ParticleFilter::resample() {

  // generate distribution according to weights
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> dist(weights.begin(), weights.end());

  // create resampled particles
  vector<Particle> resampled_particles;
  resampled_particles.resize(num_particles);

  // resample the particles according to weights
  for(int i=0; i<num_particles; i++){
    int idx = dist(gen);
    resampled_particles[i] = particles[idx];
  }

  // assign the resampled_particles to the previous particles
  particles = resampled_particles;

  // clear the weight vector for the next round
  weights.clear();
}


 

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}