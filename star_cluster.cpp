
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <pthread.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#ifdef SSE
#include <x86intrin.h>
#endif

pthread_mutex_t step_mutex;
//pthread_mutex_t sync_mutex;
pthread_barrier_t barr;

const int N          = 16;
const unsigned int MaxThreads = 4;
const long int SimSteps  = 10000000;
const int reduction  = 50;
const int dim        = 3;
float dt            = 1.0e-5;
float init_dt       = dt;
float dt_square     = 0.0f;
float accum_time    = 0;
float max_ratio     = 0; 
const float hf      = 0.5;
unsigned int interframe_counter = 0;

unsigned int lock_step;
unsigned int all_finished = (1 << MaxThreads) - 1;

typedef struct{
 unsigned int t_num;
} T_Data;

typedef struct{
  float pos[N][dim];
  float vel[N][dim];
  float acc[N][dim];
  float F[N][dim];
} Coords;

typedef struct{
  float pos[N][dim];
  float vel[N][dim];
} PosVelOnly;

typedef struct{
  float KE;
  float PE;
  float tot_E;
  float dtt;
  float KE_A[MaxThreads];
  float PE_A[MaxThreads];
  float tot_E_A[MaxThreads];
} EnergyStr;


Coords* const       T1 = new Coords;
Coords* const       T2 = new Coords;
std::vector<EnergyStr> EnergyVec;
std::vector<PosVelOnly> PosVelOnlyVec;

void MoveToCM(Coords& currPos){
  /*Make sure that the position and velocity of the center of mass of the cluster are set to zero, 
   otherwise it will be very difficult to visualize the cluster as it is zipping along. */
  float posCM[dim];
  float velCM[dim];
  for(int i = 0; i < dim; ++i){
    posCM[i] = 0;
    velCM[i] = 0;
    for(int j = 0; j < N; ++j){
      posCM[i] += currPos.pos[j][i];
      velCM[i] += currPos.vel[j][i];
    }
    posCM[i] =  posCM[i]/N;
    velCM[i] =  velCM[i]/N;
  }
  for(int i = 0; i < dim; ++i){
    for(int j = 0; j < N; ++j){
      currPos.pos[j][i] = currPos.pos[j][i] - posCM[i];
      currPos.vel[j][i] = currPos.vel[j][i] - velCM[i];
      currPos.acc[j][i] = 0;
      currPos.F[j][i]   = 0;
    }
  }
}

void adaptive_step(Coords& currPos){
  float abs_vel[N];
  float min_dist[N];
  float tmp_dist, max_dist_vel_ratio, min_dist_vel_ratio;
  int min_index;

#ifdef SSE
  /************************** SSE code ***************************/
  float distance_buffer[4] __attribute__ ((aligned (16)));
#endif

  for(int j = 0; j < N; ++j){
    abs_vel[j] = 0;  
    for(int k = 0; k < dim; ++k){
      abs_vel[j] += pow(currPos.pos[j][k], 2);
    }
    abs_vel[j] = sqrt(abs_vel[j]);

    min_dist[j] = 0;    
    
#ifdef SSE
    /************************** SSE code ***************************/
    distance_buffer[0] = currPos.pos[j][0] - currPos.pos[(j+1)%N][0];
    distance_buffer[1] = currPos.pos[j][1] - currPos.pos[(j+1)%N][1];
    distance_buffer[2] = currPos.pos[j][2] - currPos.pos[(j+1)%N][2];
    distance_buffer[3] = 0.0f;
    
    __m128 v_1 = _mm_load_ps(distance_buffer);
    
    __m128 prod = _mm_mul_ps(v_1, v_1);
    //Now sum the results
    prod = _mm_hadd_ps(prod, prod);
    prod = _mm_hadd_ps(prod, prod);
    __m128 sqrt_dist = _mm_sqrt_ss(prod);
    _mm_store_ps(distance_buffer, sqrt_dist);
    min_dist[j] = distance_buffer[0];
#else
    /************************** Non-SSE Code ***********************/
    for(int k = 0; k < dim; ++k){
      min_dist[j] += pow(currPos.pos[j][k] - currPos.pos[(j+1)%N][k], 2); //Take the distance to the next particle circularly
    }
    min_dist[j] = sqrt(min_dist[j]);
    /***************************************************************/
#endif

    for(int i = 0; i < N; ++i){
      if(i != j){
	tmp_dist = 0;
	for(int k = 0; k < dim; ++k){
	  tmp_dist += pow(currPos.pos[j][k] - currPos.pos[i][k], 2);
	}
	tmp_dist = sqrt(tmp_dist);
	if(tmp_dist < min_dist[j]){
	  min_dist[j] = tmp_dist;
	}
      }
    }    
  }
  
  max_dist_vel_ratio = min_dist[0]/abs_vel[0];
  min_dist_vel_ratio = min_dist[0]/abs_vel[0];
  min_index = 0;
  for(int j = 1; j < N; ++j){
    if( max_dist_vel_ratio < (min_dist[j]/abs_vel[j]) ){
      max_dist_vel_ratio = (min_dist[j]/abs_vel[j]);
      //min_index = j;
    }
    if( min_dist_vel_ratio > (min_dist[j]/abs_vel[j]) ){
      min_dist_vel_ratio = (min_dist[j]/abs_vel[j]);
      min_index = j;
    }
  }
  if((min_dist_vel_ratio*1.0e-4) < dt){
    dt = dt/2.0;
  }
  else if(dt < (min_dist_vel_ratio*1.0e-5)){
    dt = dt*2.0;
  } 
  max_ratio = max_dist_vel_ratio;
}

void evolve_Verlet1(T_Data& currThread){
  /************** Compute x(t + \delta t) and v(t + \delta t/2) ****************/
  float dist = 0;
  int lLim   = static_cast<int>((N*currThread.t_num)/MaxThreads);
  int uLim   = static_cast<int>((N*(currThread.t_num+1))/MaxThreads);

  dt_square  = dt*dt;

#ifdef SSE
  /************************** SSE code ***************************/
  float distance_buffer[4] __attribute__ ((aligned (16)));
  float position_buffer[4] __attribute__ ((aligned (16)));
  float velocity_buffer[4] __attribute__ ((aligned (16)));
  float accel_buffer[4]    __attribute__ ((aligned (16)));
#endif

  for(int i = lLim; i < uLim; ++i){//Process only one piece of the data allocated to the thread
    dist = 0;
    for(int k = 0; k < dim; ++k){
      T1->F[i][k] = 0;
    }
    for(int j = 0; j < N; ++j){//Need to consider interaction with all other particles
      if(i != j){
	dist = 0;

#ifdef SSE
	/************************** SSE code ***************************/
	//unroll loop
	distance_buffer[0] = T1->pos[i][0] - T1->pos[j][0];
	distance_buffer[1] = T1->pos[i][1] - T1->pos[j][1];
	distance_buffer[2] = T1->pos[i][2] - T1->pos[j][2];
	distance_buffer[3] = 0.0f;

	__m128 v_1 = _mm_load_ps(distance_buffer);
	__m128 prod = _mm_mul_ps(v_1, v_1);
	prod = _mm_hadd_ps(prod, prod);
	prod = _mm_hadd_ps(prod, prod);

	__m128 inv_sqrt_dist = _mm_rsqrt_ss(prod);
	__m128 inv_dist      = _mm_mul_ps(inv_sqrt_dist, inv_sqrt_dist);
	__m128 inv_dist_cube = _mm_mul_ps(inv_dist     , inv_sqrt_dist);
	
	_mm_store_ps(distance_buffer, inv_dist_cube);
#else
	/************************** Non-SSE Code ***********************/
	for(int k = 0; k < dim; ++k){
	  dist += pow(T1->pos[i][k] - T1->pos[j][k], 2);
	}
	dist = sqrt(dist);
	/***************************************************************/
#endif
	for(int k = 0; k < dim; ++k){
#ifdef SSE
	  T1->F[i][k] += - (T1->pos[i][k] - T1->pos[j][k])*distance_buffer[0];
#else
	  T1->F[i][k] += - (T1->pos[i][k] - T1->pos[j][k])/pow(dist, 3);
#endif	  
	}
      }
    }
#ifdef SSE
    /************************** SSE code ***************************/
    //unroll the loop
    T1->acc[i][0] = T1->F[i][0];
    T1->acc[i][1] = T1->F[i][1];
    T1->acc[i][2] = T1->F[i][2];
    
    position_buffer[0] = T1->pos[i][0];
    position_buffer[1] = T1->pos[i][1];
    position_buffer[2] = T1->pos[i][2];
    position_buffer[3] = 0.0f;

    velocity_buffer[0] = T1->vel[i][0];
    velocity_buffer[1] = T1->vel[i][1];
    velocity_buffer[2] = T1->vel[i][2];
    velocity_buffer[3] = 0.0f;

    accel_buffer[0] = T1->acc[i][0];
    accel_buffer[1] = T1->acc[i][1];
    accel_buffer[2] = T1->acc[i][2];
    accel_buffer[3] = 0.0f;

    __m128 v_pos    = _mm_load_ps(position_buffer);
    __m128 v_vel    = _mm_load_ps(velocity_buffer);
    __m128 v_acc    = _mm_load_ps(accel_buffer);
    __m128 v_dt     = _mm_load1_ps(&dt);
    __m128 v_dt_sq  = _mm_load1_ps(&dt_square);
    __m128 v_hf     = _mm_load1_ps(&hf);
    
    __m128 v_temp_pos2 = _mm_add_ps(v_pos, _mm_add_ps(_mm_mul_ps(v_vel, v_dt), _mm_mul_ps(v_hf, _mm_mul_ps(v_acc, v_dt_sq))));
    _mm_store_ps(position_buffer, v_temp_pos2);
    T2->pos[i][0] = position_buffer[0];
    T2->pos[i][1] = position_buffer[1];
    T2->pos[i][2] = position_buffer[2];
    
    __m128 v_temp_vel1 = _mm_add_ps(v_vel, _mm_mul_ps(v_hf , _mm_mul_ps(v_acc , v_dt)));
    _mm_store_ps(velocity_buffer, v_temp_vel1);
    T1->vel[i][0] = velocity_buffer[0];
    T1->vel[i][1] = velocity_buffer[1];
    T1->vel[i][2] = velocity_buffer[2];
#else
    /************************** Non-SSE Code ***********************/
    for(int k = 0; k < dim; ++k){
      T1->acc[i][k] = T1->F[i][k];
      T2->pos[i][k] = T1->pos[i][k] + T1->vel[i][k]*dt + 0.5*T1->acc[i][k]*dt_square;
      T1->vel[i][k] +=  0.5*T1->acc[i][k]*dt; //This will be v(t + \delta t/2)
    }
    /***************************************************************/
#endif
  }  
}

void evolve_Verlet2(T_Data& currThread){
  /************** Compute now v(t + \delta t) ****************/
  float dist = 0;
  int lLim   = static_cast<int>((N*currThread.t_num)/MaxThreads);
  int uLim   = static_cast<int>((N*(currThread.t_num+1))/MaxThreads);

  dt_square  = dt*dt;

#ifdef SSE
  /************************** SSE code ***************************/
  float distance_buffer[4] __attribute__ ((aligned (16)));
  float velocity_buffer[4] __attribute__ ((aligned (16)));
  float accel_buffer[4]    __attribute__ ((aligned (16)));
#endif

  for(int i = lLim; i < uLim; ++i){//Process only one piece of the data allocated to the thread
    dist = 0;
    for(int k = 0; k < dim; ++k){
      T2->F[i][k] = 0;
    }
    for(int j = 0; j < N; ++j){//Need to consider interaction with all other particles
      if(i != j){
	dist = 0;
#ifdef SSE
	/************************** SSE code ***************************/
	//unroll loop	
	distance_buffer[0] = T2->pos[i][0] - T2->pos[j][0];
	distance_buffer[1] = T2->pos[i][1] - T2->pos[j][1];
	distance_buffer[2] = T2->pos[i][2] - T2->pos[j][2];
	distance_buffer[3] = 0.0f;

	__m128 v_1 = _mm_load_ps(distance_buffer);
	__m128 prod = _mm_mul_ps(v_1, v_1);

	prod = _mm_hadd_ps(prod, prod);
	prod = _mm_hadd_ps(prod, prod);

	__m128 inv_sqrt_dist = _mm_rsqrt_ss(prod);
	__m128 inv_dist      = _mm_mul_ps(inv_sqrt_dist, inv_sqrt_dist);
	__m128 inv_dist_cube = _mm_mul_ps(inv_dist     , inv_sqrt_dist);
	
	_mm_store_ps(distance_buffer, inv_dist_cube);
#else	
	/************************** Non-SSE Code ***********************/
	for(int k = 0; k < dim; ++k){
	  dist += pow(T2->pos[i][k] - T2->pos[j][k], 2);
	}
	dist = sqrt(dist);
	/***************************************************************/
#endif
	for(int k = 0; k < dim; ++k){
#ifdef SSE	  
	  T2->F[i][k] += - (T2->pos[i][k] - T2->pos[j][k])*distance_buffer[0];
#else		  
	  T2->F[i][k] += - (T2->pos[i][k] - T2->pos[j][k])/pow(dist, 3);
#endif	  
	}
      }
    }
#ifdef SSE
    /************************** SSE code ***************************/
    //unroll the loop
    T2->acc[i][0] = T2->F[i][0];
    T2->acc[i][1] = T2->F[i][1];
    T2->acc[i][2] = T2->F[i][2];
    
    velocity_buffer[0] = T1->vel[i][0];
    velocity_buffer[1] = T1->vel[i][1];
    velocity_buffer[2] = T1->vel[i][2];
    velocity_buffer[3] = 0.0f;

    accel_buffer[0] = T2->acc[i][0];
    accel_buffer[1] = T2->acc[i][1];
    accel_buffer[2] = T2->acc[i][2];
    accel_buffer[3] = 0.0f;

    __m128 v_vel    = _mm_load_ps(velocity_buffer);
    __m128 v_acc    = _mm_load_ps(accel_buffer);
    __m128 v_dt     = _mm_load1_ps(&dt);
    __m128 v_hf     = _mm_load1_ps(&hf);
    
    __m128 v_temp_vel2 = _mm_add_ps(v_vel, _mm_mul_ps(v_hf , _mm_mul_ps(v_acc , v_dt)));
    _mm_store_ps(velocity_buffer, v_temp_vel2);
    T2->vel[i][0] = velocity_buffer[0];
    T2->vel[i][1] = velocity_buffer[1];
    T2->vel[i][2] = velocity_buffer[2];
#else
    /************************** Non-SSE Code ***********************/
    for(int k = 0; k < dim; ++k){
      T2->acc[i][k] = T2->F[i][k];
      T2->vel[i][k] = T1->vel[i][k] + 0.5*T2->acc[i][k]*dt;
    }
    /***************************************************************/
#endif
  }
}

void SumThreadEnergy(T_Data& currThread,  Coords* const T,  EnergyStr& th_E){
  float local_KE = 0;
  float local_PE = 0;
  float local_E  = 0;
  float dist = 0;

  int lLim   = static_cast<int>((N*currThread.t_num)/MaxThreads);
  int uLim   = static_cast<int>((N*(currThread.t_num+1))/MaxThreads);

#ifdef SSE
  /************************** SSE code ***************************/
  float distance_buffer[4] __attribute__ ((aligned (16)));
  float velocity_buffer[4] __attribute__ ((aligned (16)));
#endif

  for(int i = lLim; i < uLim; ++i){//Process only one piece of the data allocated to the thread

#ifdef SSE
    /************************** SSE code ***************************/
    velocity_buffer[0] = T->vel[i][0];
    velocity_buffer[1] = T->vel[i][1];
    velocity_buffer[2] = T->vel[i][2];
    velocity_buffer[3] = 0.0f;

    __m128 v_vel       = _mm_load_ps(velocity_buffer);
    __m128 v_hf        = _mm_load1_ps(&hf);
    __m128 v_temp_KE   = _mm_mul_ps(v_hf , _mm_mul_ps(v_vel , v_vel));

    v_temp_KE          = _mm_hadd_ps(v_temp_KE, v_temp_KE);
    v_temp_KE          = _mm_hadd_ps(v_temp_KE, v_temp_KE);

    _mm_store_ps(velocity_buffer, v_temp_KE);
    local_KE += velocity_buffer[0];
#else
    /************************ Non-SSE code *************************/
    for(int k = 0; k < dim; ++k){
      local_KE += 0.5*pow(T->vel[i][k], 2);
    }
    /***************************************************************/
#endif

    for(int j = i+1; j < N; ++j){//Need to consider interaction with all other particles      
      dist = 0;
#ifdef SSE
      /************************** SSE code ***************************/
      //unroll loop	
      distance_buffer[0] = T->pos[i][0] - T->pos[j][0];
      distance_buffer[1] = T->pos[i][1] - T->pos[j][1];
      distance_buffer[2] = T->pos[i][2] - T->pos[j][2];
      distance_buffer[3] = 0.0f;
      
      __m128 v_1 = _mm_load_ps(distance_buffer);
      __m128 prod = _mm_mul_ps(v_1, v_1);
     
      prod = _mm_hadd_ps(prod, prod);
      prod = _mm_hadd_ps(prod, prod);
      
      __m128 inv_dist = _mm_rsqrt_ss(prod);
      
      _mm_store_ps(distance_buffer, inv_dist);
      local_PE += -distance_buffer[0];   
#else
      /************************ Non-SSE code *************************/
      for(int k = 0; k < dim; ++k){
	dist += pow(T->pos[i][k] - T->pos[j][k], 2);
      }
      dist = sqrt(dist);
      local_PE += -1.0/dist;     
      /***************************************************************/
#endif
    }
  }
  local_E = local_PE + local_KE;

  th_E.PE_A[currThread.t_num]    = local_PE;
  th_E.KE_A[currThread.t_num]    = local_KE;
  th_E.tot_E_A[currThread.t_num] = local_E;  
}

void SumTotEnergy(EnergyStr& E_Total){
  
  float local_KE = 0;
  float local_PE = 0;
  float local_E  = 0;
  float dist = 0;
  
#ifdef E_PER_THREAD
  for(unsigned int q = 0; q < MaxThreads; ++q){
    local_PE += E_Total.PE_A[q];
    local_KE += E_Total.KE_A[q];
    local_E  += E_Total.tot_E_A[q];
  }
#endif

  for(int i = 0; i < N; ++i){//Process only one piece of the data allocated to the thread
    for(int k = 0; k < dim; ++k){
      local_KE += 0.5*pow(T1->vel[i][k], 2);
    }
    for(int j = i+1; j < N; ++j){//Need to consider interaction with all other particles      
      dist = 0;
      for(int k = 0; k < dim; ++k){
	dist += pow(T1->pos[i][k] - T1->pos[j][k], 2);
	}
      dist = sqrt(dist);
      local_PE += -1.0/dist;      
    }
  }
  local_E = local_PE + local_KE;
  
  E_Total.PE    = local_PE;
  E_Total.KE    = local_KE;
  E_Total.tot_E = local_E;
  E_Total.dtt   = dt;
}

void syncCoords(T_Data& currThread, Coords* const Ta, Coords* const Tb){
  int lLim   = static_cast<int>((N*currThread.t_num)/MaxThreads);
  int uLim   = static_cast<int>((N*(currThread.t_num+1))/MaxThreads);
  
  //pthread_mutex_lock(&sync_mutex);
  for(int si = lLim; si < uLim; ++si){//Process only one piece of the data allocated to the thread
    for(int sk = 0; sk < dim; ++sk){
      Ta->pos[si][sk] = Tb->pos[si][sk];
      Ta->vel[si][sk] = Tb->vel[si][sk];
    }
  }
  //pthread_mutex_unlock(&sync_mutex);
}

void copyCoords(PosVelOnly& Ta, Coords* const Tb){
  
  //pthread_mutex_lock(&sync_mutex);
  for(int si = 0; si < N; ++si){//Copy all data, only after all threads have finished.
    for(int sk = 0; sk < dim; ++sk){
      Ta.pos[si][sk] = Tb->pos[si][sk];
      Ta.vel[si][sk] = Tb->vel[si][sk];
    }
  }
  //pthread_mutex_unlock(&sync_mutex);
}

void WriteOutCoords(std::string f_out, std::vector<PosVelOnly>& CVec){
  static int junk_1 = 0;
  std::ofstream fout;
  fout.open(f_out.c_str());
  if (fout.is_open()){
    fout.precision(3);
    for(std::vector<PosVelOnly>::iterator C_Pos = CVec.begin(); C_Pos != CVec.end(); ++C_Pos){
      fout << junk_1 << " ";
      for(int j = 0; j < N; ++j){
	for(int k = 0; k < dim; ++k){
	  fout << C_Pos->pos[j][k] << " ";
	}
	for(int k = 0; k < dim; ++k){
	  fout << C_Pos->vel[j][k] << " ";
	}	
      }        
      fout << std::endl;
      fout.flush();
      ++junk_1;
    }
  }
  fout.close();
  CVec.clear();
}

void WriteOutEnergy(std::string f_out, std::vector<EnergyStr>& EVec){  
  int counter = 0;

  std::ofstream fout;
  fout.open(f_out.c_str());
  if (fout.is_open()){
    for(std::vector<EnergyStr>::iterator C_E = EVec.begin(); C_E != EVec.end(); ++C_E){
      fout << counter    << " "  
	   << C_E->KE    << " " 
	   << C_E->PE    << " " 
	   << C_E->tot_E << " " 
	   << C_E->dtt   << " "
	   << max_ratio  << " "
	   << std::endl;
      ++counter;
    }
  }
  fout.close();
  EVec.clear();
}

void Tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " "){
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos){
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
  }
}

void readInitCoords(std::string f_in, Coords* const T){
  
  std::ifstream fin;
  std::string c_line;
  std::vector<std::string> line_elements;
  int line = 0;
  int col  = 0;
  fin.open(f_in.c_str());
  if (fin.is_open()){
    while( fin ){
      getline(fin, c_line);
      Tokenize(c_line, line_elements);      
      
      col = 0;
      for(std::vector<std::string>::iterator v_elem = line_elements.begin(); v_elem != line_elements.end(); ++v_elem){
	if(col < dim){
	  T->pos[line][col] = atof(v_elem->c_str());
	}
	else{
	  T->vel[line][col%dim] = atof(v_elem->c_str());
	}
	++col;
      }
      line_elements.clear();
      ++line;
    }
  }
  fin.close();
}

void* F(void* argg){
  T_Data tmp_argg( *( static_cast<T_Data*> (argg) ) );
  pthread_mutex_lock(&step_mutex);
  lock_step = lock_step | (1 << tmp_argg.t_num);
  pthread_mutex_unlock(&step_mutex);

  for(int k = 0; k < SimSteps; ++k){      
    /**************** Do the real work here ********************/

    if(k%3 == 0){
      evolve_Verlet1(tmp_argg);
    }
    else if(k%3 == 1){
      evolve_Verlet2(tmp_argg);
    }
    else if(k%3 == 2){
      syncCoords(tmp_argg, T1, T2);
    }
    
    /***********************************************************/
    pthread_mutex_lock(&step_mutex);
    // Aquire lock and check if all other threads have finished.
    lock_step = lock_step | (1 << tmp_argg.t_num);        

    if(lock_step != all_finished){//Not the last thread, wait for the others to finish.   
      //Unlock the step mutex before waiting, otherwise the lock is never released.
      pthread_mutex_unlock(&step_mutex);
      int rc = pthread_barrier_wait(&barr);
      if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD){
	std::cout << "Could not wait on barrier" << std::endl;
	exit(-1);
      }
    }
    else{
      /* All other threads have finished and are waiting.
	 We can now compute things like energy, changes in time step, 
	 and check if enough time has accumulated to write out the coordinates*/
      lock_step = 0;    
      if( k%3 == 2 ){
	++interframe_counter;
	if(reduction <= interframe_counter){
	  interframe_counter = 0;
	  adaptive_step(*T2);
	  
	  if( 3*reduction*init_dt < accum_time){
	    PosVelOnly c_Pos;
	    EnergyStr c_E_Tot;

	    accum_time = 0;
#ifdef CALC_ENERGY
	    SumTotEnergy(c_E_Tot);
	    EnergyVec.push_back(c_E_Tot);
#endif
	    copyCoords(c_Pos, T2);
	    PosVelOnlyVec.push_back(c_Pos);
	  }
	}
	else{
	  accum_time += dt;
	}
      }
               
      //Unlock the step mutex before waiting otherwise will create a deadlock.
      pthread_mutex_unlock(&step_mutex); 
      int rc = pthread_barrier_wait(&barr);
      if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD){
	std::cout << "Could not wait on barrier" << std::endl;
	exit(-1);
      }
    }
  }
  pthread_exit(NULL);
}

int main(int argc, char* argv[]){
  pthread_attr_t attr;
  unsigned int thr_id[MaxThreads];
  T_Data* TT[MaxThreads];
  pthread_t  p_thread[MaxThreads];

  for(unsigned int i = 0; i < MaxThreads; ++i){
    TT[i] = new T_Data;
    TT[i]->t_num = i;
  }

  lock_step = 0;  
  readInitCoords("initial_data.dat", T1);
  MoveToCM(*T1);
  *T2 = *T1;
  //readInitCoords("data_d.dat", T2);
  //MoveToCM(*T2);
   
  pthread_mutex_init(&step_mutex, NULL);  
  //pthread_mutex_init(&sync_mutex, NULL);
   
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  if(pthread_barrier_init(&barr, NULL, MaxThreads)){
    std::cout << "Could not create a barrier" << std::endl;
    exit(-1);
  }

  /* Allow any thread to run on any CPU */
#ifdef LINUX
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  for (unsigned int j = 0; j < MaxThreads; j++){
    CPU_SET(j, &cpuset);
  }
#endif
  for(unsigned int i = 0; i < MaxThreads; ++i){
    /* Force each thread to run on only one CPU. */
#ifdef LINUX
    CPU_ZERO(&cpuset);
    CPU_SET(i, &cpuset);
#endif
    thr_id[i] = pthread_create(&(p_thread)[i], &attr, F, static_cast<void*>(TT[i]));
#ifdef LINUX
    int s = pthread_setaffinity_np(p_thread[i], sizeof(cpu_set_t), &cpuset);
    if (s != 0){
      std::cerr << "pthread_setaffinity_np() , s = " << s << std::endl;
    }

    s = pthread_getaffinity_np(p_thread[i], sizeof(cpu_set_t), &cpuset);
    printf("Set returned by pthread_getaffinity_np() contained:\n");
    for (int j = 0; j < CPU_SETSIZE; j++){
      if (CPU_ISSET(j, &cpuset)){
	printf(" CPU %d\n", j);
      }
    }
#endif
  }
  
  for(unsigned int i = 0; i < MaxThreads; ++i){
    pthread_join(p_thread[i], NULL);
  }
  std::cout << "Threads finished " << std::endl;
  fflush(stdout);

#ifdef CALC_ENERGY
  WriteOutEnergy("energy.txt", EnergyVec);
#endif
  WriteOutCoords("positions.txt", PosVelOnlyVec);

  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&step_mutex);

  //pthread_mutex_destroy(&sync_mutex);
  pthread_barrier_destroy(&barr);

  for(unsigned int i = 0; i < MaxThreads; ++i){
    delete TT[i];
  }

  pthread_exit(NULL);
}
