#ifndef _IMU_H_
#define _IMU_H_

#include "stm32f10x.h"
#include <math.h>
#include <stdint.h>

typedef struct Vel_Pos  
{
	float vx;
	float vy;
	float vz;
	float px;
	float py;
	float pz;
}VelPosPackage;

typedef struct IMU
{
	float gx;
	float gy;
	float gz;
	float ax;
	float ay;
	float az;
}IMU_filtered;

//----------------------------------------------------------------------------------------------------
// Variable declaration
extern float q0, q1, q2, q3;	// quaternion of sensor frame relative to auxiliary frame
extern float roll, pitch, yaw;
extern volatile VelPosPackage send_vel_pos;
extern IMU_filtered imu_gyro_acc;
//---------------------------------------------------------------------------------------------------
void AttitudeUpdate(void);
void QuarternionUpdate(IMU_filtered gyro_acc);
void Kalman_filter(int16_t *input);
void InerNav(float ax,float ay,float az);
void cal_vel_pos(float *dacc);

#endif
