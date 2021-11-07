#include <iostream>
#include <memory>
#include <cmath>
#include <vector>
#include "../../eigen/Eigen/Dense"
#include "../../eigen/Eigen/Sparse"
#include <SFML/Graphics.hpp>

template<typename T>
void DiscreteLaplacian1D(size_t n,Eigen::SparseMatrix<T> &mat){
	std::vector<Eigen::Triplet<T>> triplets(4+(n-2)*3);
	int k=0;
	triplets[0] = Eigen::Triplet<T>(0,0,-1.);
	++k;
	triplets[k] = Eigen::Triplet<T>(0,1,2);
	++k;
	for(int i = 1; i < n-1; ++i){
		triplets[k] = Eigen::Triplet<T>(i,i-1,-1);
		++k;
		triplets[k] = Eigen::Triplet<T>(i,i,2);
		++k;
		triplets[k] = Eigen::Triplet<T>(i,i+1,-1);
		++k;
	}	
	
	
	triplets[k] = Eigen::Triplet<T>(n-1,n-2,-1);
	++k;	
	triplets[k] = Eigen::Triplet<T>(n-1,n-1,2);
	mat.setFromTriplets(triplets.begin(),triplets.end());
}

template<typename T>
void DiscreteLaplacian2D(size_t n,Eigen::SparseMatrix<T> &mat){
	std::vector<Eigen::Triplet<T>> triplets;
	for(int i = 0; i < n*n; ++i){
		int j = i;
		triplets.push_back(Eigen::Triplet<T>(i,j,4));
//		mat.insert(i,j) = 4;		
		
		j = i+1;
		if(j <= n*n - 1)  triplets.push_back(Eigen::Triplet<T>(i,j,-1.));
		
			
		j = i+n;
		if(j <= n*n - 1) triplets.push_back(Eigen::Triplet<T>(i,j,-1.));
		
		
		j = i-1;
		if(j >= 0) triplets.push_back(Eigen::Triplet<T>(i,j,-1.));
		
			
		j = i-n;
		if(j >= 0) triplets.push_back(Eigen::Triplet<T>(i,j,-1.));
	}

	mat.setFromTriplets(triplets.begin(), triplets.end());
}

void euler(Eigen::VectorXf &u,const Eigen::VectorXf &F, float dt){
		u = u + F*dt;

}	


template<typename T>
void F_calc(Eigen::VectorXf &F,const Eigen::SparseMatrix<T> &L,const Eigen::VectorXf &u, float dt, float t){
	F = L*u;
	t = t+dt;	

}


void get_image(unsigned char * buffer, const Eigen::VectorXf data, const size_t n){
	for(size_t i=0; i<n; ++i){
		unsigned char val = (unsigned char)(255*(1+data[i])*0.5);
		buffer[i*4]   = val;	
		buffer[i*4+1] = val;
		buffer[i*4+2] = 255;
		buffer[i*4+3] = 255;
	}

}

void get_starting_point(Eigen::VectorXf& u, const size_t n){
	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			u(i*n + j) = sin(0.05*(double)i) *  sin(0.05*(double)j);
		}
	}
}
int main(int argc, char** argv)
{
	size_t n{400};

	for(int i=0; i < argc; ++i){ std::cout << argv[i] << " ";}
	std::cout << std::endl;
	if(argc > 1) n = atoi(argv[1]);
	std::cout << argc << " " << n << std::endl;
	size_t width{n};
	size_t height{n};	
	float t{0}, dt{0.001};

	//allocating solution vector to calculate at each frame
	Eigen::VectorXf u(n*n);
	Eigen::VectorXf F(n*n);

	//allocating a matrix w x h and building the laplacian
	Eigen::SparseMatrix<float> L(width*height,width*height);
	L.reserve(5*width*height);
	DiscreteLaplacian2D(n,L);

	//setting up sfml anti aliasing settings, better quality but worse performance
	sf::ContextSettings settings;
	settings.antialiasingLevel = 8;

   	sf::RenderWindow window(sf::VideoMode(width,height),
		       	"Heat simulation",
			sf::Style::Default,settings);
	

	//setting up the buffer to build the image	
    	unsigned char *buffer = new unsigned char[width*height*4];
//	std::cout << u.cols() << "" << u.rows() << std::endl;
	get_starting_point(u,n);
	get_image(buffer,u,n*n);


	//create the image to display
	sf::Texture texture;
	texture.create(width,height);
    
    	sf::Sprite sprite(texture);
	texture.update(buffer);
	window.display();
	bool once{1};	
	while (window.isOpen())
    	{
        sf::Event event;

        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::KeyPressed){
                if(once){
                    once = 0;
                    printf("Key pressed \n");
		    get_starting_point(u,n);
                    texture.update(buffer);
                }

            }
            if (event.type == sf::Event::Closed){
                window.close();
            }
            if(event.type == sf::Event::KeyReleased){
                printf("Key released\n");
                window.clear();
            }
            if(event.type == sf::Event::Resized){
                window.clear();
               // width = window.getSize().x;
               // height = window.getSize().y;
                texture.update(buffer);
            }
       	 }
	for(int i = 0; i < 5; ++i){
		F_calc(F,L,u,dt,t);
		euler(u,F,dt);
	}
	get_image(buffer,u,n*n);	
	texture.update(buffer);
        window.draw(sprite);

        window.display();

    }
}
