#include "imu.h"
#include "mpu6050.h"
#include "usart.h"
#include <math.h>

//---------------------------------------------------------------------------------------------------
// Variable definitions
#define EarthRadius  6371393.0f
#define EarthOmega  0.00007292f
#define gravity 9.8015f    //北京地区重力加速度
#define pi 3.1415926f
#define gyroconst 0.00106422f //(pi/180)*Data/16.4 (ras/s)
//#define accconst 0.00119629f  //A=9.8*Data/8192  (m/s^2)
#define accconst 0.00012207f  //A=Data/8192  (g)
#define dt 0.005f //采样周期
#define Kp 0.85f                        
#define Ki 0.002f                        
#define halfT 0.0025f    //采样周期的一半           
float exInt = 0.0, eyInt = 0.0, ezInt = 0.0;   
float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;	// quaternion of sensor frame relative to auxiliary frame
float roll = 0.0f, pitch = 0.0f, yaw = 0.0f; 
float T11=0.0f,T12=0.0f,T13=0.0f;
float T21=0.0f,T22=0.0f,T23=0.0f;
float T31=0.0f,T32=0.0f,T33=0.0f;

float gyro_offset[3]={7.063431f,43.826538f,-45.191451f};
float acc_offset[3]={0.052172f,-0.027099f,-0.131154f};//离线matlab球体拟和求得
float acc_sens[3]={1.008307f,1.010849f,0.981117f};

#define Q 0.005f
#define R 0.75f 

float Kk[3]={0.0,0.0,0.0};
float Pk[3]={0.0,0.0,0.0};
float Pk_1[3]={0.0,0.0,0.0};
float Xk[3]={0.0,0.0,0.0};
float Xk_1[3]={0.0,0.0,0.0};
volatile VelPosPackage send_vel_pos;   
IMU_filtered imu_gyro_acc;

#define SMOOTHLEN 30L
#define INDEX   6L
int16_t buffer[INDEX][SMOOTHLEN]={0};

int16_t smooth(int16_t data,uint8_t index)
{
	uint8_t i;
	long long int sum=0;
	int16_t average;
	
	for(i=SMOOTHLEN-1;i>0;i--)
	{
		buffer[index][i]=buffer[index][i-1];
	}
	buffer[index][0]=data;
	
	for(i=0;i<SMOOTHLEN;i++)
	{
		sum+=buffer[index][i];
	}
	average=(int16_t)(sum/SMOOTHLEN);
	
	return average;
	
}
float gyro_w1[3]={0.0},gyro_w2[3]={0.0};

void AttitudeUpdate(void)  
{ 
	  uint8_t i;
    int16_t acc[3],gyr[3],acc_smooth[3],gyro_smooth[3];
	  float buf_gyro[3],buf_acc[3];
    
	  mpu6050_read_accel(&acc[0], &acc[1], &acc[2]);
    mpu6050_read_gyro(&gyr[0], &gyr[1], &gyr[2]);

		for(i=0;i<3;i++)
		{
			acc_smooth[i]=smooth(acc[i],i);
      gyro_smooth[i]=smooth(gyr[i],i+3);
		}
		for(i=0;i<3;i++)
		{
			gyro_w1[i]=(float)gyro_smooth[i]*0.3f;
      gyro_w2[i]= gyro_w2[i]*0.7f+gyro_w1[i];
		}
    Kalman_filter(acc_smooth);
	
	  for(i=0;i<3;i++)
    {
				buf_gyro[i]=(gyro_w2[i]-gyro_offset[i])*gyroconst; //rad/s
				buf_acc[i]=(Xk_1[i]*accconst-acc_offset[i])*acc_sens[i];  //g 
    }
		
		imu_gyro_acc.gx=buf_gyro[0];
		imu_gyro_acc.gy=buf_gyro[1];
		imu_gyro_acc.gz=buf_gyro[2];
		imu_gyro_acc.ax=buf_acc[0];
		imu_gyro_acc.ay=buf_acc[1];
		imu_gyro_acc.az=buf_acc[2];
		
		QuarternionUpdate(imu_gyro_acc); 
}
                                
void QuarternionUpdate(IMU_filtered gyro_acc) 
{
  float q0temp,q1temp,q2temp,q3temp;
	float gx, gy, gz;
	float ax, ay, az;
  float vx, vy, vz;
  float ex, ey, ez;      
	float norm; 
	
  float q0q0 = q0*q0;
  float q0q1 = q0*q1;
  float q0q2 = q0*q2;
	float q0q3 = q0*q3;
	
  float q1q1 = q1*q1;
	float q1q2 = q1*q2;
  float q1q3 = q1*q3;
	
  float q2q2 = q2*q2;
  float q2q3 = q2*q3;
  float q3q3 = q3*q3;  
	
	
	gx=gyro_acc.gx;
	gy=gyro_acc.gy;
	gz=gyro_acc.gz;
	
	ax=gyro_acc.ax;
	ay=gyro_acc.ay;
	az=gyro_acc.az;
	
	//计算坐标变换矩阵Cb2n
	T11=1.0f-2.0f*(q2q2+q3q3);
	T12=2.0f*(q1q2-q0q3);
	T13=2.0f*(q1q3+q0q2);

	T21=2.0f*(q1q2+q0q3);
	T22=1.0f-2.0f*(q1q1+q3q3);
	T23=2.0f*(q2q3-q0q1);

	T31=2.0f*(q1q3-q0q2);
	T32=2.0f*(q2q3+q0q1);
	T33=1.0f-2.0f*(q1q1+q2q2);
	
  if(ax*ay*az==0)return;
	
  norm = sqrt(ax*ax + ay*ay + az*az);
  ax = ax /norm;
  ay = ay / norm;
  az = az / norm;
  //将重力加速度[0,0,1]'投影到载体坐标系
  vx = T31;        
  vy = T32;
  vz = T33;
  //误差计算
  ex = (ay*vz - az*vy) ;                                                                  
  ey = (az*vx - ax*vz) ;
  ez = (ax*vy - ay*vx) ;
	//计算误差积分
  exInt = exInt + ex * Ki;                                           
  eyInt = eyInt + ey * Ki;
  ezInt = ezInt + ez * Ki;
  // 比例积分控制器
  gx = gx + Kp*ex + exInt;  
  gy = gy + Kp*ey + eyInt;
  gz = gz + Kp*ez + ezInt;    
  
  q0temp=q0;
  q1temp=q1;
  q2temp=q2;
  q3temp=q3;
  //一阶毕卡法
  q0 = q0temp + (-q1temp*gx - q2temp*gy -q3temp*gz)*halfT;
  q1 = q1temp + (q0temp*gx + q2temp*gz -q3temp*gy)*halfT;
  q2 = q2temp + (q0temp*gy - q1temp*gz +q3temp*gx)*halfT;
  q3 = q3temp + (q0temp*gz + q1temp*gy -q2temp*gx)*halfT;
  //四元数规范化处理
  norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
  q0 = q0 / norm;
  q1 = q1 / norm;
  q2 = q2 / norm;
  q3 = q3 / norm;
	
	//yaw+=-gz*dt*57.3f;
	//pitch = asin(T32)*57.3f; // pitch
  //roll = atan2(-T31,T33)*57.3f; // roll
}


void Kalman_filter(int16_t *input)
{ 
	int i=0;
	for(i=0;i<3;i++)
	{
		Xk[i] = Xk_1[i];
		Pk[i] =Pk_1[i] + Q;
		Kk[i] = Pk[i]/(Pk[i] + R); 
		Xk_1[i] =Xk[i] + Kk[i]*((float)input[i]-Xk[i]);
		Pk_1[i] = (1.0f -Kk[i]) * Pk[i];
	}
}

struct Vector3f
{
	float x;	
	float y;
	float z;
};


struct Vector3f  Cb2n(struct Vector3f data)
{
	struct Vector3f out;
	
	out.x=T11*data.x+T12*data.y+T13*data.z;
	out.y=T21*data.x+T22*data.y+T23*data.z;
	out.z=T31*data.x+T32*data.y+T33*data.z;
	
	return out;
}


struct Vector3f  Cn2b(struct Vector3f data)
{
	struct Vector3f out;
	
	out.x=T11*data.x+T21*data.y+T31*data.z;
	out.y=T12*data.x+T22*data.y+T32*data.z;
	out.z=T13*data.x+T23*data.y+T33*data.z;
	
	return out;
}


struct Vector3f crossproduct(struct Vector3f V1,struct Vector3f V2)
{
	struct Vector3f out;
	
	out.x=V1.y*V2.z-V1.z*V2.y;
	out.y=V1.z*V2.x-V1.x*V2.z;
	out.z=V1.x*V2.y-V1.y*V2.x;
	
	return out;
}

#define Kpos 0.01f
#define Kvel 0.1f
float  pos[3]={0.0f,0.0f,0.0f};
float  vel[3]={0.0f,0.0f,0.0f};

void cal_vel_pos(float *dacc)
{
		
		if(dacc[0]>-1.0f&&dacc[0]<1.0f)dacc[0]=0.0f;
		if(dacc[1]>-1.0f&&dacc[1]<1.0f)dacc[1]=0.0f;
		if(dacc[2]>-1.0f&&dacc[2]<1.0f)dacc[2]=0.0f;
	
	
	  //速度积分得位移
		vel[0]+=Kvel*dacc[0]*dt;
		vel[1]+=Kvel*dacc[1]*dt;
		vel[2]+=Kvel*dacc[2]*dt;

		pos[0]+=Kpos*vel[0]*dt;
		pos[1]+=Kpos*vel[1]*dt;
		pos[2]+=Kpos*vel[2]*dt;
}


const struct Vector3f grav={0.0f,0.0f,1.0f};
#define DEEP 10L
float buf_facc[3][DEEP]={0.0f};
double init_sum[3]={0.0f,0.0f,0.0f};
float init_offset[3]={0.0f,0.0f,0.0f};
uint16_t wait_time=0;

void InerNav(float ax,float ay,float az)
{
	uint8_t i;
	struct Vector3f bgra;
	struct Vector3f facc;
	double buf_sum[3]={0.0f,0.0f,0.0f};
	float buf_ave[3]={0.0f,0.0f,0.0f};
	float dfacc[3]={0.0f,0.0f,0.0f};
	//float dvel[3]={0.0f,0.0f,0.0f};
	
	//bgra=Cn2b(grav); 
	
	bgra.x=T31;
	bgra.y=T32;
	bgra.z=T33; 
	
	facc.x=ax-bgra.x; 
	facc.y=ay-bgra.y;
	facc.z=az-bgra.z;
	
	wait_time++;
	if(wait_time>=1000&&wait_time<=1100)
	{
		init_sum[0]+=facc.x;
		init_sum[1]+=facc.y;
		init_sum[2]+=facc.z;
	}
	if(wait_time==1101)
	{
		init_offset[0]=init_sum[0]/100.0f;
		init_offset[1]=init_sum[1]/100.0f;
		init_offset[2]=init_sum[2]/100.0f;
		for(i=0;i<DEEP;i++)
		{
			buf_facc[0][i]=init_offset[0];
			buf_facc[1][i]=init_offset[1];
			buf_facc[2][i]=init_offset[2];
		}
	}
	if(wait_time>=1102)
	{
		wait_time=1102;
		
		for(i=DEEP-1;i>0;i--)
		{
			buf_facc[0][i]=buf_facc[0][i-1];
			buf_facc[1][i]=buf_facc[1][i-1];
			buf_facc[2][i]=buf_facc[2][i-1];
		}
		buf_facc[0][0]=facc.x;
		buf_facc[1][0]=facc.y;
		buf_facc[2][0]=facc.z;
		
		for(i=0;i<DEEP;i++)
		{
			buf_sum[0]+=buf_facc[0][i];
			buf_sum[1]+=buf_facc[1][i];
			buf_sum[2]+=buf_facc[2][i];
		}
		
		buf_ave[0]=(float)(buf_sum[0]/DEEP);
		buf_ave[1]=(float)(buf_sum[1]/DEEP);
		buf_ave[2]=(float)(buf_sum[2]/DEEP);
		
		
		dfacc[0]=(buf_ave[0]-init_offset[0])*gravity;
		dfacc[1]=(buf_ave[1]-init_offset[1])*gravity;
		dfacc[2]=(buf_ave[2]-init_offset[2])*gravity;
		
		cal_vel_pos(dfacc); 
		
		//数据发送
		send_vel_pos.vx=vel[0];
		send_vel_pos.vy=vel[1];
		send_vel_pos.vz=vel[2];
		send_vel_pos.px=pos[0];
		send_vel_pos.py=pos[1];
		send_vel_pos.pz=pos[2];
		
	}
}



