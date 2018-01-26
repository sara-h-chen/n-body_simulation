// Translate this file with
//
// g++ -O3 --std=c++11 spaceboddies.c -o spacebodies
//
// Run it with
//
// ./spacebodies
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2017 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <cstdio>

struct Body {
  double mass;
  double velocityX, velocityY, velocityZ;
  double positionX, positionY, positionZ;
  double forceX, forceY, forceZ;
  bool isActive;
};

double t = 0;
double tFinal = 0;

int NumberOfBodies = 0;

Body *bodies;

std::ofstream videoFile;


// ----------------------------------------------------
//                    SETUP FUNCTION
// ----------------------------------------------------

void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-2) / 7;

  bodies = new Body[NumberOfBodies];

  int readArgument = 1;

  tFinal = std::stof(argv[readArgument]); readArgument++;

  // Process each body separately
  for (int i=0; i < NumberOfBodies; ++i) {
    // The first three numbers passed into the command line are the x variables
    // representing the position of the particle
    bodies[i].positionX = std::stof(argv[readArgument]); readArgument++;
    bodies[i].positionY = std::stof(argv[readArgument]); readArgument++;
    bodies[i].positionZ = std::stof(argv[readArgument]); readArgument++;

    // The next three numbers passed in represent the velocity vectors
    bodies[i].velocityX = std::stof(argv[readArgument]); readArgument++;
    bodies[i].velocityY = std::stof(argv[readArgument]); readArgument++;
    bodies[i].velocityZ = std::stof(argv[readArgument]); readArgument++;

    bodies[i].mass = std::stof(argv[readArgument]); readArgument++;
    bodies[i].isActive = true;

    if (bodies[i].mass<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
}


// ----------------------------------------------------
//              PARAVIEW RELATED FUNCTIONS
// ----------------------------------------------------

void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}


void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}

/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */
void printParaviewSnapshot(int counter) {
  std::stringstream filename;
  filename << "paraview/result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i < NumberOfBodies; ++i) {
    out << bodies[i].positionX
        << " "
        << bodies[i].positionY
        << " "
        << bodies[i].positionZ
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}


// ----------------------------------------------------
//                  UPDATE FUNCTIONS
// ----------------------------------------------------

void updatePosition(const double timeStepSize) {
  double accelerationX, accelerationY, accelerationZ;
  double deltaT = (timeStepSize * timeStepSize) / 2;

  for (int i=0; i < NumberOfBodies; ++i) {
    accelerationX = bodies[i].forceX / bodies[i].mass;
    accelerationY = bodies[i].forceY / bodies[i].mass;
    accelerationZ = bodies[i].forceZ / bodies[i].mass;

    bodies[i].positionX = bodies[i].positionX + (timeStepSize * bodies[i].velocityX) + (deltaT * accelerationX);
    bodies[i].positionY = bodies[i].positionY + (timeStepSize * bodies[i].velocityY) + (deltaT * accelerationY);
    bodies[i].positionZ = bodies[i].positionZ + (timeStepSize * bodies[i].velocityZ) + (deltaT * accelerationZ);

    bodies[i].velocityX = bodies[i].velocityX + (timeStepSize * accelerationX);
    bodies[i].velocityY = bodies[i].velocityY + (timeStepSize * accelerationY);
    bodies[i].velocityZ = bodies[i].velocityZ + (timeStepSize * accelerationZ);

    // Set all force to 0 for calculation of next iteration
    bodies[i].forceX = 0;
    bodies[i].forceY = 0;
    bodies[i].forceZ = 0;
  }
}

// Part 2: If you have two bodies headed towards each other then they must collide
// Part 3: Make the time step change according to how close the bodies are to one another so that the particles don't just pass through each other
void updateBody() {

  const double timeStepSize = 0.0000001;
  
  // Step 1.1: All bodies interact and move
  for (int i=0; i < NumberOfBodies; ++i) {
    for (int j=0; j < NumberOfBodies; ++j) {
      if (i != j) {
        // Distance between body i and every other body
    		const double distance = sqrt(
      		(bodies[i].positionX-bodies[j].positionX) * (bodies[i].positionX-bodies[j].positionX) +
      		(bodies[i].positionY-bodies[j].positionY) * (bodies[i].positionY-bodies[j].positionY) +
      		(bodies[i].positionZ-bodies[j].positionZ) * (bodies[i].positionZ-bodies[j].positionZ)
    		);

        // DEBUG
        // printf("\nDistance : %5.7f", distance); 
        
        // TODO: Decrease time step size as min distance gets closer 
        // TODO: If distance < 1e-8 then fuse the bodies
        if (distance < 1e-8) {
          printf("Should collide with distance : %5.2f", distance);


        }

        // } else {
          double massDistance = bodies[i].mass * bodies[j].mass / (distance * distance * distance);

          // Calculate the force between body and the others
          bodies[i].forceX += (bodies[j].positionX-bodies[i].positionX) * massDistance ;
          bodies[i].forceY += (bodies[j].positionY-bodies[i].positionY) * massDistance ;
          bodies[i].forceZ += (bodies[j].positionZ-bodies[i].positionZ) * massDistance ;
        // }

        // Increase time step
        t += timeStepSize;
      }
    }

    printf("\nBody %d: %7.8f  %7.8f  %7.8f", i, bodies[i].positionX, bodies[i].positionY, bodies[i].positionZ);
  }

  updatePosition(timeStepSize);
}


// ----------------------------------------------------
//                     MAIN METHOD
// ----------------------------------------------------

int main(int argc, char** argv) {
  // Insufficient args
  if (argc==1) {
    std::cerr << "please add the final time plus a list of object configurations as tuples px py pz vx vy vz m" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "100.0   0 0 0   1.0 0 0   1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "100.0   0 0 0   1.0 0 0   1.0     0 1.0 0   1.0 0 0   1.0 \t One spiralling around the other one" << std::endl
              << "100.0   3.0 0 0   0 1.0 0   0.4     0 0 0   0 0 0   0.2     2.0 0 0   0 0 0   1.0 \t Three body setup from first lecture" << std::endl
              << std::endl;
    return -1;
  }
  // Mismatched args
  else if ( (argc-2)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    return -2;
  }

  setUp(argc,argv);
  openParaviewVideoFile();
  printParaviewSnapshot(0);

  int timeStepsSinceLastPlot = 0;
  const int plotEveryKthStep = 100;
  while (t<=tFinal) {
    updateBody();
    timeStepsSinceLastPlot++;
    if (timeStepsSinceLastPlot%plotEveryKthStep==0) {
      // printParaviewSnapshot(timeStepsSinceLastPlot/plotEveryKthStep);
    }
  }

  // closeParaviewVideoFile();

  return 0;
}
