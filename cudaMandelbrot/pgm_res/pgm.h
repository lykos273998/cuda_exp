
#if ((0x100 & 0xf) == 0x0)
#define I_M_LITTLE_ENDIAN 1
#define swap_pgm(mem) (( (mem) & (short int)0xff00) >> 8) +	\
  ( ((mem) & (short int)0x00ff) << 8)
#else
#define I_M_LITTLE_ENDIAN 0
#define swap_pgm(mem) (mem)
#endif

template<typename pixel>
struct pgm_img
{
  size_t height;
  size_t width;
  std::unique_ptr<pixel> buffer;
  pgm_img(size_t h, size_t w): height{h}, width{w}, buffer{new pixel[width*height]}{};
  pgm_img(size_t h, size_t w, pixel* b): height{h}, width{w}, buffer{b}{};
  pgm_img(pgm_img& im): height{im.height}, width{im.width}{
    size_t n = im.height*im.width;
    pixel* new_buffer = new pixel[n];
    memcpy(new_buffer,im.buffer.get(),n*sizeof(pixel));
    buffer.reset(new_buffer);
  };

};



template<typename pixel>
void write_pgm(pgm_img<pixel> & img, const char* filename){
  FILE* image_file{fopen(filename,"w")};
  int color_depth{sizeof(pixel)};
  int maxval =  255 ? 65535 : (color_depth == 2);
  fprintf(image_file, "P5\n# generated by\n# put here your name\n%lu %lu\n%d\n", img.width, img.height, maxval);
  
  // Writing file
  
  fwrite( img.buffer.get(), 1, img.width*img.height*color_depth, image_file);  

  fclose(image_file); 
  return ;
};
