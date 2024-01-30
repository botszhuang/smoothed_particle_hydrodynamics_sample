#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct 
{
    double x, y ;
} xyStruct;

typedef struct
{ 
  xyStruct x ;
  xyStruct v ;
  double density ;
  double pressure ;
  xyStruct pressure_force  ;
  xyStruct viscosity_force ;
  xyStruct force ;
} testParticleProfileStruct ;

#define MAX_NUMBER_PARTICLES 125
#define DOMAIN_WALL   40
#define DOMAIN_HEIGHT 80
#define PARTICLE_MASS 1
#define ISOTROPIC_EXPONENT 20
#define BASE_DENSITY 1
#define SMOOTHING_LENTH  5
#define DYNAMIC_VISCOSITY 0.5
#define DAMPING_COEFFICIENT -0.9
#define TIME_STEP_LENGTH 0.01
#define N_TIME_STEPS 2000
#define ADD_PARTICELS_EVERY 50
#define PLOT_EVERY 5

const double DOMAIN_X_MIN = SMOOTHING_LENTH ;
const double DOMAIN_X_MAX = DOMAIN_WALL - SMOOTHING_LENTH ;

const double DOMAIN_Y_MIN = SMOOTHING_LENTH ;
const double DOMAIN_Y_MAX = DOMAIN_HEIGHT - SMOOTHING_LENTH ;
const double NORMALIZATION_DENSITY = 315 * PARTICLE_MASS / (  64 * 3.14159 * pow ( SMOOTHING_LENTH , 9.0 ) ) ;
const double NORMALIZATION_PRESSURE_FORCE = - 45 * PARTICLE_MASS / (  3.14159 * pow ( SMOOTHING_LENTH , 6.0 ) ) ;
const double NORMALIZATION_VISCOSITY_FORCE =  45 * DYNAMIC_VISCOSITY * PARTICLE_MASS / (  3.14159 * pow ( SMOOTHING_LENTH , 6.0 ) ) ;

const xyStruct CONSTANT_FORCE = ( xyStruct ) { .x = 0 , .y = -0.1 } ;

void initParticles ( testParticleProfileStruct * p ) ;
void printParticles ( const unsigned int currentParticleNumber ,
                      testParticleProfileStruct * p ) ;
double get_distance ( xyStruct a , xyStruct b ) ;
void calculate_distance ( const unsigned int currentParticleNumber , 
                          testParticleProfileStruct * particles  ) ;
void calculate_density  ( const unsigned int currentParticleNumber ,
                          testParticleProfileStruct * particles  ) ;
void calculate_pressure ( const unsigned int currentParticleNumber , 
                          testParticleProfileStruct * particles ) ;
void calculate_pressure_forece  (  const unsigned int currentParticleNumber ,
                                   testParticleProfileStruct * particles  ) ;
void calculate_viscosity_forece  ( const unsigned int currentParticleNumber ,
                                   testParticleProfileStruct * particles  ) ;  
void calculate_forece  (  const unsigned int currentParticleNumber ,
                          testParticleProfileStruct * particles  ) ; 
void Euler_step  (  const unsigned int currentParticleNumber ,
                    testParticleProfileStruct * particles  ) ;
void enforce_boudary_condition  (  const unsigned int currentParticleNumber ,
                                   testParticleProfileStruct * particles  ) ; 
double distanceMatrix [ MAX_NUMBER_PARTICLES ] [ MAX_NUMBER_PARTICLES ] ;

int main()
{
 printf ( "Hello!\n" ) ;
 system ( "rm particle.txt");

 printf ( "NORMALIZATION_DENSITY         = % .3e\n", NORMALIZATION_DENSITY );
 printf ( "NORMALIZATION_PRESSURE_FORCE  = % .3e\n", NORMALIZATION_PRESSURE_FORCE  );
 printf ( "NORMALIZATION_VISCOSITY_FORCE = % .3e\n", NORMALIZATION_VISCOSITY_FORCE );

 testParticleProfileStruct particles [ MAX_NUMBER_PARTICLES ] ;
 initParticles  ( particles ) ;
 printParticles ( MAX_NUMBER_PARTICLES ,  particles ) ;

 unsigned int currentParticleNumber = 0 ;
 for ( unsigned int l = 0 ; l < N_TIME_STEPS ; l++ ) {

   if ( currentParticleNumber < MAX_NUMBER_PARTICLES )
   {
    if ( ( l % ADD_PARTICELS_EVERY ) == 0 ) 
    {
        currentParticleNumber += 3 ;
        if ( currentParticleNumber > MAX_NUMBER_PARTICLES ) 
        {
        currentParticleNumber = MAX_NUMBER_PARTICLES ;
        }
    }
   }
   
   calculate_distance ( currentParticleNumber ,particles ) ;
   calculate_density  ( currentParticleNumber ,particles ) ;
   calculate_pressure ( currentParticleNumber ,particles ) ;
   calculate_pressure_forece ( currentParticleNumber ,particles ) ;
   calculate_viscosity_forece ( currentParticleNumber ,particles ) ;
   calculate_forece ( currentParticleNumber , particles ) ;

   Euler_step ( currentParticleNumber , particles ) ;
   enforce_boudary_condition ( currentParticleNumber , particles ) ; 
   
   if ( ( l % PLOT_EVERY ) == 0 ){
     printParticles ( currentParticleNumber ,  particles ) ;
   }
 }
 
 puts("making a movie ...") ;
 system ("gnuplot plot.gu");

 return EXIT_SUCCESS ;
}
void enforce_boudary_condition  (  const unsigned int currentParticleNumber ,
                                   testParticleProfileStruct * particles  ) 
{
  #define X particles[i].x
  #define V particles[i].v  
   for ( unsigned int i = 0 ; i < currentParticleNumber ; i++  ) 
   { 
    if ( X.x < DOMAIN_X_MIN ) { X.x = DOMAIN_X_MIN ; V.x *= DAMPING_COEFFICIENT ; }
    if ( X.x > DOMAIN_X_MAX ) { X.x = DOMAIN_X_MAX ; V.x *= DAMPING_COEFFICIENT ; }
    if ( X.y < DOMAIN_Y_MIN ) { X.y = DOMAIN_Y_MIN ; V.y *= DAMPING_COEFFICIENT ; }
    if ( X.y > DOMAIN_Y_MAX ) { X.y = DOMAIN_Y_MAX ; V.y *= DAMPING_COEFFICIENT ; }
   }
  #undef X
  #undef V
}
void Euler_step (  const unsigned int currentParticleNumber ,
                   testParticleProfileStruct * particles  ) 
{
  #define X particles[i].x
  #define V particles[i].v  
  #define force   particles[i].force 
  #define density   particles[i].density
   for ( unsigned int i = 0 ; i < currentParticleNumber ; i++  ) 
   { 
    V.x += TIME_STEP_LENGTH * force.x / density ;
    V.y += TIME_STEP_LENGTH * force.y / density ; 
    X.x += TIME_STEP_LENGTH * V.x ;
    X.y += TIME_STEP_LENGTH * V.y ;                   
   }
  #undef X
  #undef V
  #undef force 
  #undef density
}                          
void calculate_forece  (  const unsigned int currentParticleNumber ,
                          testParticleProfileStruct * particles  ) 
{
   for ( unsigned int i = 0 ; i < currentParticleNumber ; i++  ) 
   { 
    particles[i].force.x = CONSTANT_FORCE.x 
                       + particles[i].pressure_force.x 
                       + particles[i].viscosity_force.x ;
    particles[i].force.y = CONSTANT_FORCE.y
                       + particles[i].pressure_force.y 
                       + particles[i].viscosity_force.y ; 
   }
}                          
void calculate_viscosity_forece  ( const unsigned int currentParticleNumber ,
                                   testParticleProfileStruct * particles  ) 
{
   testParticleProfileStruct host , guest ;
   double distance ;
   double Ld1 ;

   for ( unsigned int i = 0 ; i < currentParticleNumber ; i++  ) 
   { 
     host = particles [i] ;
     host.viscosity_force = ( xyStruct ) { 0.0 , 0.0 } ;
     for ( unsigned int j = 0 ; j < currentParticleNumber ; j++  ) 
     {  
        if ( i == j ) { continue ; }
        distance = distanceMatrix [i][j] ;
        if ( distance > SMOOTHING_LENTH ) { continue ; }
        guest = particles [j] ;

        Ld1  = SMOOTHING_LENTH - distance ;
        Ld1 /= guest.density ;

        host.viscosity_force.x += ( guest.v.x - host.v.x ) * Ld1  ;
        host.viscosity_force.y += ( guest.v.y - host.v.y ) * Ld1  ;
     }
     particles [i].viscosity_force.x = NORMALIZATION_VISCOSITY_FORCE * host.viscosity_force.x ;
     particles [i].viscosity_force.y = NORMALIZATION_VISCOSITY_FORCE * host.viscosity_force.y ;
    }

}                           
void calculate_pressure_forece  ( const unsigned int currentParticleNumber ,
                                  testParticleProfileStruct * particles  ) 
{
   testParticleProfileStruct host , guest ;
   double distance ;
   double Ld2 ;

   for ( unsigned int i = 0 ; i < currentParticleNumber ; i++  ) 
   { 
     host = particles [i] ;
     host.pressure_force = ( xyStruct ) { 0.0 , 0.0 } ;
     for ( unsigned int j = 0 ; j < currentParticleNumber ; j++  ) 
     {  
        if ( i == j ) { continue ; }
        distance = distanceMatrix [i][j] ;
        if ( distance > SMOOTHING_LENTH ) { continue ; }
        guest = particles [j] ;

        Ld2  = SMOOTHING_LENTH - distance ;
        Ld2 *= Ld2 ;
        Ld2  = Ld2 / distance
                * ( guest.pressure + host.pressure ) 
                / ( 2 * guest.density ) ;

        host.pressure_force.x -= ( guest.x.x - host.x.x ) * Ld2  ;
        host.pressure_force.y -= ( guest.x.y - host.x.y ) * Ld2  ;
     }
     particles [i].pressure_force.x = NORMALIZATION_PRESSURE_FORCE * host.pressure_force.x ;
     particles [i].pressure_force.y = NORMALIZATION_PRESSURE_FORCE * host.pressure_force.y ;
    }
}
void calculate_pressure ( const unsigned int currentParticleNumber , testParticleProfileStruct * particles )
{
   for ( unsigned int i = 0 ; i < currentParticleNumber ; i++  ) 
   { 
      particles [i].pressure = ISOTROPIC_EXPONENT * ( particles[i].density - BASE_DENSITY ) ;
   }
}
void calculate_density  ( const unsigned int currentParticleNumber , testParticleProfileStruct * particles  ) 
{
  double Ld3 ;
  double distance ;
  double density ;
   for ( unsigned int i = 0 ; i < currentParticleNumber ; i++  ) 
   { 
      density = 0 ;
      Ld3 = 0 ;
      for ( unsigned int j = 0 ; j < currentParticleNumber ; j++ ) 
      {
        distance = distanceMatrix [i][j] ;
        if ( distance > SMOOTHING_LENTH ) { continue ; }

        Ld3 = ( ( SMOOTHING_LENTH * SMOOTHING_LENTH ) - ( distance * distance ) ) ; 
        density += ( Ld3 * Ld3 * Ld3 ) ;          
      }
      particles [i].density = NORMALIZATION_DENSITY * density ;
   }
}
void calculate_distance ( const unsigned int currentParticleNumber , testParticleProfileStruct * particles  ) 
{
   for ( unsigned int i = 0 ; i < currentParticleNumber ; i++  ) {
   for ( unsigned int j = 0 ; j < currentParticleNumber ; j++  ) {
        if ( i == j ) { continue; }
        distanceMatrix[i][j] = get_distance ( particles [i].x , particles [j].x ) ;
   }}
}
void printParticles ( const unsigned int currentParticleNumber ,
                      testParticleProfileStruct * p ) 
{
    static unsigned int piccounter = 0 ;

    FILE *fp ;
    char * fname = "particle.txt" ;
    
    fp = fopen( fname , "a");
    if (!fp) {
        printf("cannot open file, %s...\n", fname ) ;
        exit(1);
    }
       

    fprintf ( fp , "#%u) X, Y , Vx, Vy\n", piccounter );
    for ( unsigned int i = 0 ; i < currentParticleNumber ; i++ ) 
    {
        fprintf ( fp , "% 10.6e, % 10.6e , % 10.6e, % 10.6e \n" , 
                  p[i].x.x , p[i].x.y , p[i].v.x , p[i].v.y );
    }
    
    fprintf ( fp , "\n\n" ) ;
    fflush ( fp ) ;
    fclose ( fp ) ;

    piccounter ++ ;

}
double get_distance ( xyStruct a , xyStruct b ) {
    return sqrt( ( ( a.x - b.x ) * ( a.x - b.x ) ) + ( ( a.y - b.y ) * ( a.y - b.y ) )  ) ;
}
void initParticles ( testParticleProfileStruct * p ) 
{
    srand ( 123 ) ;
    const double x0 = DOMAIN_WALL / 8.0 ;
    double x = 0 ; 
    double r = 0 ;
    for ( unsigned int i = 0 ; i < MAX_NUMBER_PARTICLES ; i++ ) 
    {
        r = (double) rand() / (RAND_MAX + 1.0);  
        if      ( ( i % 3 ) == 0 ) { x = DOMAIN_X_MIN + ( x0 ) ; }
        else if ( ( i % 3 ) == 1 ) { x = DOMAIN_X_MIN + ( x0 +   SMOOTHING_LENTH      ) ; }
        else                       { x = DOMAIN_X_MIN + ( x0 + ( SMOOTHING_LENTH *2 ) ) ; }
        p[i].x = ( xyStruct ) { x + r , DOMAIN_Y_MAX } ;
        p[i].v = ( xyStruct ) { -3.0 , -15.0 } ;
        p[i].density = 0 ;
        p[i].pressure = 0 ;
    }
}