#include <stdio.h>
#include <iostream>
#include <complex>
#include <memory>
#include <cstdlib>
#include <cmath>
#include <pgm.h>
#include <ctime>
#include <omp.h>

#include <fstream>
#include <cpu_info.hpp>
#include <cuda_runtime.h>


#define NITER 1000
#define THRESHOLD 1000

#define MAXVAL 65535




int powi (int base, unsigned int exp)
{
    int res = 1;
    while (exp) {
        if (exp & 1)
            res *= base;
        exp >>= 1;
        base *= base;
    }
    return res;
}

__device__
void mset_calc(unsigned short &k, double i, double j, double scale, double ofx, double ofy){
  /**
   * Device code used to perform the actual calculation of each pixel
   * */
    double Re, Im;
    Re = i*scale + ofx;
    Im = j*scale + ofy;
    double zRe{0.}, zIm{0.}; 
    double cRe{Re},cIm{Im}, z_1Re, z_1Im, temp_zRe, temp_zIm;

    double t;
    z_1Re = zRe;
    z_1Im = zIm;


    zRe = cRe;
    zIm = cIm;
    k = 1;
    while(k < NITER){
        //
        temp_zRe = zRe;
        temp_zIm = zIm;

        t = zRe*zRe - zIm*zIm + z_1Re + cRe;
        zIm = 2*zRe*zIm + z_1Im + cIm;

        zRe = t;
        z_1Re = temp_zRe;
        z_1Im = temp_zIm;
        if(zRe*zRe + zIm*zIm > THRESHOLD){
            break;
        }
        ++k;
        
  }


}

void HOST_mset_calc(unsigned short &k, double i, double j, double scale, double ofx, double ofy){
    double Re, Im;
    Re = (i + ofx)*scale;
    Im = (j + ofy)*scale;
    double zRe{0.}, zIm{0.}; 
    double cRe{Re},cIm{Im}, z_1Re, z_1Im, temp_zRe, temp_zIm;

    double t;
    z_1Re = zRe;
    z_1Im = zIm;


    zRe = cRe;
    zIm = cIm;
    k = 1;
    while(k < NITER){
        //
        temp_zRe = zRe;
        temp_zIm = zIm;

        t = zRe*zRe - zIm*zIm + z_1Re + cRe;
        zIm = 2*zRe*zIm + z_1Im + cIm;

        zRe = t;
        z_1Re = temp_zRe;
        z_1Im = temp_zIm;
        if(zRe*zRe + zIm*zIm > THRESHOLD){
            break;
        }
        ++k;
        
  }


}


__global__
void CudaMandelbrot(size_t width, size_t height,unsigned short* image, double scale, double cx, double cy)
{
  /**
   * Kernel used to act like a bridge between CUDA and CPU code
   * The strange thing about cuda is that, every "action" on the loop
   * is performed by a single CUDA_thread, so each CUDA_thread executes a single
   * simple task. The GPU seems to be more efficient when under full load
  */
  size_t i{blockIdx.x*blockDim.x + threadIdx.x};
  size_t j{blockIdx.y*blockDim.y + threadIdx.y};

  if(i < width && j < height){
    double ii{ (double)(i)/height - 0.5};
    double jj{ (double)(j)/width - 0.5};
    unsigned short k{0};
    mset_calc(k,ii,jj,scale,cx,cy);
    image[i*width + j] = k;
  }
}

void print_err_msg(cudaError & err){
  if (err != cudaSuccess)
    {
        fprintf(stdout, "%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char** argv)
{
  //int N{255 ? 65535 : 1};

  //default hw
  size_t height{1200};
  size_t width{1200};
  clock_t t;

  double scale{0.01}, cx{-0.7}, cy{-0.5};
  if(argc > 4){
    //argv for position of the fractal
    scale = atof(argv[2]);
    cx = atof(argv[3]);
    cy = atof(argv[4]);
  }

  if(argc > 5){
    //argv for height and width
    height=atoll(argv[5]);
    width=atoll(argv[6]);
  }

  //allocate buffer on host
  unsigned short* myimg = new unsigned short[height*width];
  std::cout << "Img dimensions: h = " << height << "  w = " << width << "\n";
  std::cout << "Paramenters: scale = "<< scale << " cx = " << cx << " cy = " << cy << std::endl;

  /**everything follows is compiled only in "profiling mode, allocates a new buffer for 
   * storing the image calculated on the cpu, not needed if you do not use that
    */
  #ifdef PROF
    unsigned short* myimg_CPU = new unsigned short[height*width];
    std::cout << "Running mandelbrot (like) set calculation with time profiling: CPU vs GPU" << std::endl;
   // std::cout << "using threads on CPU " << get_cpu_info() << std::endl;
    
    t = clock();
    #pragma omp parallel
    {
      #pragma omp single
      {
        std::cout << "Using " << omp_get_num_threads() << " threads on CPU " << get_cpu_info() << std::endl;
      }
      #pragma omp for
      for(int i=0; i< height; ++i){
        for(int j=0; j< width; ++j){
          double ii{ (double)(i)/height - 0.5};
          double jj{ (double)(j)/width - 0.5};
          unsigned short k{0};
          HOST_mset_calc(k,ii,jj,scale,cx,cy);
          myimg_CPU[i*width + j] = k;
        }
      }
    }
    std::cout << "*** Elapsed calculation time CPU: " << (double)(clock() - t)/CLOCKS_PER_SEC << std::endl;
    std::cout  << "Printing a number only to trick gcc to compile" << " " << myimg_CPU[1] << std::endl;

  #endif


  
  cudaError err = cudaSuccess;

  //generating number of threads to spawn on GPU
  dim3 threadsPerBlock(16, 16);
  dim3 numBlocks(width / threadsPerBlock.x, height / threadsPerBlock.y);
  
  //allocate array on GPU and check for errors
  unsigned short* gpuImg;
  err = cudaMalloc((void**)&gpuImg,width*height*2);

  print_err_msg(err);

  //let GPU calculate mandelbrot set and check for errors
  
  t = clock();

  CudaMandelbrot<<<numBlocks, threadsPerBlock>>>(width, height, gpuImg, scale, cx, cy);

  
  
  err = cudaGetLastError();
  std::cout << "*** Elapsed calculation time GPU: " << (double)(clock() - t)/CLOCKS_PER_SEC << std::endl;

  

  print_err_msg(err);

  //copy back pgm image
  t = clock();
  err = cudaMemcpy(myimg,gpuImg,height*width*2, cudaMemcpyDeviceToHost);
  std::cout << "*** Time to copy the result from GPU: " << (double)(clock() - t)/CLOCKS_PER_SEC << std::endl;
  print_err_msg(err);

  //check for img being equal
  /*
  #ifdef PROF
  for(int i = 0; i< height*width;++i){
    if(myimg[i] != myimg_CPU[i]){
      std::cout << "err" << myimg[i] - myimg_CPU[i] << std::endl;
    }
  }
  #endif
  */
  //wirte pgm image on file
  pgm_img<unsigned short> test_img{height,width,myimg};
  
  std::cout << "writing image on file ";

  if(argc > 1){
    write_pgm(test_img, argv[1]);
    std::cout << argv[1] << std::endl;
    
  }
  else{
    write_pgm(test_img, "img.ppm");
    std::cout << "img.pgm" << std::endl;
  }

  cudaFree(gpuImg);
  std::cout << "End" << std::endl;
  
  return 0;
}