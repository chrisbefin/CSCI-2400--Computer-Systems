#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <fstream>
#include "Filter.h"
#include <omp.h>

using namespace std;

#include "rdtsc.h"

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);

int
main(int argc, char **argv)
{

  if ( argc < 2) {
    fprintf(stderr,"Usage: %s filter inputfile1 inputfile2 .... \n", argv[0]);
  }

  //
  // Convert to C++ strings to simplify manipulation
  //
  string filtername = argv[1];

  //
  // remove any ".filter" in the filtername
  //
  string filterOutputName = filtername;
  string::size_type loc = filterOutputName.find(".filter");
  if (loc != string::npos) {
    //
    // Remove the ".filter" name, which should occur on all the provided filters
    //
    filterOutputName = filtername.substr(0, loc);
  }

  Filter *filter = readFilter(filtername);

  double sum = 0.0;
  int samples = 0;

  for (int inNum = 2; inNum < argc; inNum++) {
    string inputFilename = argv[inNum];
    string outputFilename = "filtered-" + filterOutputName + "-" + inputFilename;
    struct cs1300bmp *input = new struct cs1300bmp;
    struct cs1300bmp *output = new struct cs1300bmp;
    int ok = cs1300bmp_readfile( (char *) inputFilename.c_str(), input);

    if ( ok ) {
      double sample = applyFilter(filter, input, output);
      sum += sample;
      samples++;
      cs1300bmp_writefile((char *) outputFilename.c_str(), output);
    }
    delete input;
    delete output;
  }
  fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

struct Filter *
readFilter(string filename)
{
  ifstream input(filename.c_str());

  if ( ! input.bad() ) {
    int size = 0;
    input >> size;
    Filter *filter = new Filter(size);
    int div;
    input >> div;
    filter -> setDivisor(div);
    for (int i=0; i < size; i++) {
      for (int j=0; j < size; j++) {
        int value;
        input >> value;
        filter -> set(i,j,value);
      }
    }
    return filter;
  } else {
    cerr << "Bad input in readFilter:" << filename << endl;
    exit(-1);
  }
}


double
applyFilter(struct Filter *filter, cs1300bmp *input, cs1300bmp *output)
{
    long long cycStart, cycStop;
    cycStart = rdtscll();

    int columns = (input->width) -  1;
    int rows = (input->height) - 1;
    output -> width = columns + 1;
    output -> height = rows + 1;
    
    int size = filter -> getSize();//code motion
    int acc1, acc2, acc3, acc4, acc5, acc6, acc7, acc8, acc9 = 0;
    float divisor = filter -> getDivisor();//code motion
    int *arrayStart = (filter -> get(0, 0));//using direct memory access to eliminate a function call in the inner loop
    int i = 0;
    
    #pragma omp parallel for
    for(int row = 1; row < rows; row++) {
        for(int col = 1; col < columns; col++) {
            acc1 = 0;
            acc2 = 0;
            acc3 = 0;
            acc4 = 0;
            acc5 = 0;
            acc6 = 0;
            acc7 = 0;
            acc8 = 0;
            acc9 = 0;
            
            acc1 += (input->color[0][row+i-1][col-1]*arrayStart[i*size]);
            acc2 += (input->color[0][row+i][col-1]*arrayStart[(i+1)*size]);
            acc3 += (input->color[0][row+i+1][col-1]*arrayStart[(i+2)*size]);
            acc4 += (input->color[1][row+i-1][col-1]*arrayStart[i*size]);
            acc5 += (input->color[1][row+i][col-1]*arrayStart[(i+1)*size]);
            acc6 += (input->color[1][row+i+1][col-1]*arrayStart[(i+2)*size]);
            acc7 += (input->color[2][row+i-1][col-1]*arrayStart[i*size]);
            acc8 += (input->color[2][row+i][col-1]*arrayStart[(i+1)*size]);
            acc9 += (input->color[2][row+i+1][col-1]*arrayStart[(i+2)*size]);
                    
            acc1 += (input->color[0][row+i-1][col]*arrayStart[i*size+1]);
            acc2 += (input->color[0][row+i][col]*arrayStart[(i+1)*size+1]);
            acc3 += (input->color[0][row+i+1][col]*arrayStart[(i+2)*size+1]);
            acc4 += (input->color[1][row+i-1][col]*arrayStart[i*size+1]);
            acc5 += (input->color[1][row+i][col]*arrayStart[(i+1)*size+1]);
            acc6 += (input->color[1][row+i+1][col]*arrayStart[(i+2)*size+1]);
            acc7 += (input->color[2][row+i-1][col]*arrayStart[i*size+1]);
            acc8 += (input->color[2][row+i][col]*arrayStart[(i+1)*size+1]);
            acc9 += (input->color[2][row+i+1][col]*arrayStart[(i+2)*size+1]);
                    
            acc1 += (input->color[0][row+i-1][col+1]*arrayStart[i*size+2]);
            acc2 += (input->color[0][row+i][col+1]*arrayStart[(i+1)*size+2]);
            acc3 += (input->color[0][row+i+1][col+1]*arrayStart[(i+2)*size+2]);
            acc4 += (input->color[1][row+i-1][col+1]*arrayStart[i*size+2]);
            acc5 += (input->color[1][row+i][col+1]*arrayStart[(i+1)*size+2]);
            acc6 += (input->color[1][row+i+1][col+1]*arrayStart[(i+2)*size+2]);
            acc7 += (input->color[2][row+i-1][col+1]*arrayStart[i*size+2]);
            acc8 += (input->color[2][row+i][col+1]*arrayStart[(i+1)*size+2]);
            acc9 += (input->color[2][row+i+1][col+1]*arrayStart[(i+2)*size+2]);             
            
            acc1 = (acc1+acc2+acc3) / divisor;
            acc4 = (acc4+acc5+acc6) / divisor;
            acc7 = (acc7+acc8+acc9) / divisor;
            
            (acc1 < 0) ? acc1=0 : acc1=acc1;
            (acc1 > 255) ? acc1=255 : acc1=acc1;
            
            (acc4 < 0) ? acc4=0 : acc4=acc4;
            (acc4 > 255) ? acc4=255 : acc4 = acc4;
            
            (acc7 < 0) ? acc7=0 : acc7=acc7;
            (acc7 > 255) ? acc7=255 : acc7=acc7;
            
            output -> color[0][row][col] = acc1;
            output -> color[1][row][col] = acc4;
            output -> color[2][row][col] = acc7;
        }
    }
  cycStop = rdtscll();
  double diff = cycStop - cycStart;
  double diffPerPixel = diff / (output -> width * output -> height);
  fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n",
	  diff, diff / (output -> width * output -> height));
  return diffPerPixel;
}