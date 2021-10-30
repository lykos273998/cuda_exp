#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<omp.h>
#include<unistd.h>
#include<math.h>

int main(int argc, char** argv){
    
    long long int npx, npy, npoints;
    unsigned int niter = 1000;
    double threshold = 1000.;
    
    if(argc < 2){
        printf("insert a number of point to calculate the mandelbrot set\n");
        printf("if you want insert also the offset of the region to calculate\n the default is Re(z) [0,1] Im(z) [0,1]");
        return 0;
    }

    FILE *out;
    out = fopen("m_set_v3.dat", "w");
    //fprintf(out, "norm\n");
    npx = atoll(argv[1]);
    npy = atoll(argv[2]);
    double cx = 0.;
    double cy = 0.;
    double scale = 1.;
    if(argc > 5){
        cx = atof(argv[3]);
        cy = atof(argv[4]);
        scale = atof(argv[5]);
    }

    npoints = npx*npy;
    long long int ppt = npoints;
    unsigned int *it_to_explode;
    it_to_explode = (unsigned int *)calloc(ppt, sizeof(unsigned int));
    #pragma omp parallel
    {
        #pragma omp master
        {
            #ifdef _OPENMP
            printf("using %d threads\nperforming %lld sampling points with center [%lf %lf] and scale %lf \n", omp_get_num_threads(), npoints, cx ,cy, scale);
            #endif
        }
        #pragma omp barrier

        #ifdef _OPENMP
            printf("Thread %d starting! \n", omp_get_thread_num());
        #endif
        
        //double *cRE, *cIM;

        
        //cRE = (double*)malloc(ppt*sizeof(double));
        //cIM = (double*)malloc(ppt*sizeof(double));
        long long int i_max = 0;
        long long int i = 0;
        double complex z, c, z_1, temp_z;
        double Re, Im;
        int stride = 16;
        i_max = stride* npoints/stride;
        long long int i_ext;


        #pragma omp for nowait schedule(dynamic)
        for(i_ext = 0; i_ext < npoints; i_ext += stride)
            for(i = i_ext; i < i_ext + stride; i++ ){
                
                Re = (i % npx)/((double)npx) - 0.5;
                Im = (i / npx)/((double)npy) - 0.5;
                Re = Re*scale + cx;
                Im = Im*scale + cy;
                c = Re + I*Im;
                z = 0. + 0.*I;

                z_1 = z;
                z = c;
                unsigned int k = 1;
                for(k = 2; k < niter; k++){
                    //
                    temp_z = z;
                    z = z*z + z_1 + c;
                    z_1 = temp_z;
                    if(creal(conj(z)*z) > threshold){
                        break;
                    }
                    
                }

            
            it_to_explode[i] = k;
            
        }

       
        //remainder
        #pragma omp for nowait schedule(dynamic)
        for(i = i_max; i < npoints; i++){
                
                Re = (i % npx)/((double)npx) - 0.5;
                Im = (i / npx)/((double)npy) - 0.5;
                Re = Re*scale + cx;
                Im = Im*scale + cy;
                c = Re + I*Im;
                z = 0. + 0.*I;

                z_1 = z;
                z = c;
                unsigned int k = 1;
                for(k = 2; k < niter; k++){
                    //
                    temp_z = z;
                    z = z*z + z_1 + c;
                    z_1 = temp_z;
                    if(creal(conj(z)*z) > threshold){
                        break;
                    }
                    
                }

            
            it_to_explode[i] = k;
            
        }
        




        
        #ifdef _OPENMP
            printf("Thread %d finished! \n", omp_get_thread_num());
        #endif
        
        
        
    }



    for(long long int i = 0; i< npy; i++){
        for(long long int j = 0; j<npx; j++){
        fprintf(out, "%u ",it_to_explode[i*npx+j]);
        }
        fprintf(out, "\n");

    }



    printf("exported to dat file!\n");
    //int o = system("python3 draw.py");
    return 0;

}