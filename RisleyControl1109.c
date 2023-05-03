//define lib

#include <stdio.h> 
#include <string.h> 
#include <sys/types.h> 
#include <errno.h> 
#include <sys/stat.h> 
#include <fcntl.h> 
#include <unistd.h> 
#include <termios.h> 
#include <stdlib.h> 
 #include<math.h>
 #include<stdlib.h>
 #include<pthread.h>
#include<time.h>

#define NUM_THREADS 6
#define n 5         //define the Akima length
int systimes = 0;
double x_00[n], y_00[n],x_01[n],y_01[n];
clock_t start, stop1,stop2,stop3;
double time_doing_1,time_doing_2;
unsigned char buff4[1];
unsigned char buff5[1];
unsigned char buff6[1];

// ******************************************************************************
//SerialPortsVaries
int fd4,fd5,fd6;
int get4,get5,get6;
int timer=0;
unsigned char TC[13];
unsigned char TC_1[50];
unsigned char TC_2[50];
unsigned char enable[3]={'e','n',0x0d};
unsigned char unable[2]={'k',0x0d};
char chaxunweizhi[4]={'p','f','b',0x0D};
//TrackVaries
unsigned char trackStatus=0x00;
//TuoBaLiang
int detaX=0,detaY=0;
char tuoBaLiangStatus=0x00;
//DianJiShiJiWeiZhi
double pos1=0;
double pos2=0;
double L0;
//KongZhiCanShu
// double kRotary=1,kAngle=1,kvv=50,KI0=0.1,kp=1,ki=0,kd=0;
double kRotary=0.65,kAngle=0.65,kvv=50,KI0=0.1,kp=1,ki=0,kd=0.1;
/*data *********************************************************************************/
double PI=3.141592653;
double ddetaXa1=0,ddetaXa2=0,ddetaXa3=0,ddetaXa4=0,ddetaXa5=0,ddetaXa6=0,ddetaXa7=0,ddetaXa8=0,ddetaXa9=0,ddetaXa0=0;
double ddetaYa0 =0,ddetaYa1 =0,ddetaYa2 =0,ddetaYa3 =0,ddetaYa4 =0,ddetaYa5 =0,ddetaYa6 =0,ddetaYa7 =0,ddetaYa8 =0,ddetaYa9 =0;
// int solutionPro_status =0,trackStatus=0,enableStatus=0,solutionStatus =0;
// double k1=0.3,k2=0.3,k3=0.3,k4=0.3,pos1=0,pos2=0;
// char *buff4[100],*buff5,*buff6;
// char buf[] = "00000110";


// double a1= 18.15 / 180 * 3.141595265358979323846264338327952884197169399;
// double  a2 = 18.15 / 180 * 3.14159265358979323846264338327952884197169399;
// double n1 = 1.551;
// double n2 = 1.551;
// double zeroPoint = 0.000000000000001;
// double limit = 0.447;
// double xFOV = 3.1826;//3.85
// double yFOV = 2.5461;//2.375
// double xPIxel = 1280;
// double yPIxel = 1024;
// double x[2];
// double forwardSolution[2];
// double inverseSolution[5];
// double targetPos[2];
// double zyTheta[2];



 //define hanshu############################################################

void move_behind(double data[] ,double n_data,int j) {
	for (int i = 0; i < j-1; i++) {
		data[i] = data[i + 1];
	}
	data[j-1] = n_data;
} 

double* Akima_forecast(double x,double y, double nx,int d) {
	double x_0[n];
	double y_0[n]; 
	double m[n+6];
	double t[n+1];
	double k[4];
	if (d == 0) {
		move_behind(x_00, x, n);
		move_behind(y_00, y, n);
		for (int i = 0; i < n; i++) {
			x_0[i] = x_00[i];
			y_0[i] = y_00[i];
		}
	}
	else {
		move_behind(x_01, x, n);
		move_behind(y_01, y, n);
		for (int i = 0; i < n; i++) {
			x_0[i] = x_01[i];
			y_0[i] = y_01[i];
		}
	}
	if (systimes < n+1) {
		static double fun[2];
		fun[0] = x;
		fun[1] = y;
		return fun;
	}
	else {
		for (int i = 0; i < n - 1; i++) {
			m[i] = (y_0[i + 1] - y_0[i]) / (x_0[i + 1] - x_0[i]);
		}
		m[n-1] = (m[n - 2] + m[n - 3]) / 2;
		m[n] = (m[n-1] + m[n - 2]) / 2;
		m[n+1] = (m[n] + m[n - 1]) / 2;
		m[n + 2] = (m[n+1] + m[n]) / 2;
		m[n + 3] = (m[n+2] + m[n+1]) / 2;
		m[n + 4] = (m[n+3] + m[n+2]) / 2;
		for (int i = 0; i < n+1; i++) {
			if (m[i] == m[i + 1] && m[i+2] == m[i + 3]) {
				t[i] = (m[i + 1] + m[i+2]) / 2;
			}
			else if (m[i + 4] == m[i + 3] && m[i+2] == m[i + 1]) {
				t[i] = (m[i+2] + m[i + 3]) / 2;
			}
			else {
				t[i] = (fabs(m[i + 3] - m[i+2]) * m[i + 1] + fabs(m[i + 1] - m[i]) * m[i+2]) / (fabs(m[i + 3] - m[i+2]) + fabs(m[i + 1] - m[i]));
			}
		}
		k[0] = y_0[n - 1];
		k[1] = t[n - 1];
		k[2] = (3 * m[n - 1] - 2 * t[n - 1] - t[n]) / (x_0[n - 1] - x_0[n - 2]);
		k[3] = (t[n] + t[n - 1] - 2 * m[n - 1]) / ((x_0[n] - x_0[n - 1])* (x_0[n] - x_0[n - 1]));
		static double fun[2];
		fun[0] = nx;
		fun[1] = k[0] + k[1] * (nx - x_0[n - 1]) + k[2] * (nx - x_0[n - 1]) * (nx - x_0[n - 1]) + k[3] * (nx - x_0[n - 1]) * (nx - x_0[n - 1]) * (nx - x_0[n - 1]);
		return fun;
	}
}

int CharArrayToInt(unsigned char char1,unsigned char char2,unsigned char char3,unsigned char char4)
{
    int ret_val=0;
    ret_val=char4;
    ret_val<<=8;
    ret_val|=char3;
    ret_val<<=8;
    ret_val|=char2;
    ret_val<<=8;
    ret_val|=char1;
    return ret_val;
}

unsigned char * IntToCharArry(int data)
{
    static unsigned char charArray[4];
    charArray[0]=(char)data;
    data=data>>8;
    charArray[1]=(char)data;
    data=data>>8;
    charArray[2]=(char)data;
    data=data>>8;
    charArray[3]=(char)data;
    data=data>>8;
    return charArray;
}

//define ShowCharArryOX
void ShowCharArryOX(char data[],int j){
    for(int i =0;i<j;i++){
        printf("%0x\n",data[i]);
    }
}

//define mod
double mod(double a,double b)
{
    double c;
    c=fmod(a,b);
    if (c<0){c=c+b;}
    return c;
}

//define convertToString
char* ToString(double dater)
{
    static char s[20];
    sprintf(s,"%.4f",dater);
    return s;
}


// define it
unsigned char ChecksumBit8(unsigned char data[], int len){
    unsigned char dataSum=0;
    for(int j =0;j<len;j++){
        dataSum +=data[j];}
        return dataSum;
}


 //define  open_port
int open_port(int fd,int comport) 
{ 
	char *dev[]={"/dev/ttyS3","/dev/ttyS4","/dev/ttyS5"}; 

	if (comport==4)
	{
		fd = open( "/dev/ttyS3", O_RDWR|O_NOCTTY|O_NDELAY); 
		if (-1 == fd)
		{ 
			perror("Can't Open Serial Port"); 
			return(-1); 
		} 
     } 
     else if(comport==5)
     {     
		fd = open( "/dev/ttyS4", O_RDWR|O_NOCTTY|O_NDELAY); 
		if (-1 == fd)
		{ 
			perror("Can't Open Serial Port"); 
			return(-1); 
		} 
     } 
     else if (comport==6)
     { 
		fd = open( "/dev/ttyS5", O_RDWR|O_NOCTTY|O_NDELAY); 
		if (-1 == fd)
		{ 
			perror("Can't Open Serial Port"); 
			return(-1); 
		} 
     } 

     if(fcntl(fd, F_SETFL, 0)<0) 
     		printf("fcntl failed!\n"); 
     else 
		printf("fcntl=%d\n",fcntl(fd, F_SETFL,0)); 

     if(isatty(STDIN_FILENO)==0) 
		printf("standard input is not a terminal device\n"); 
     else 
		printf("isatty success!\n"); 
        printf("fd-open=%d\n",fd); 
     return fd; 
}

//define  set_port

int set_opt(int fd,int nSpeed, int nBits, char nEvent, int nStop) 
{ 
     struct termios newtio,oldtio; 

     if  ( tcgetattr( fd,&oldtio)  !=  0) {  
      perror("SetupSerial 1");
	printf("tcgetattr( fd,&oldtio) -> %d\n",tcgetattr( fd,&oldtio)); 
      return -1; 
     } 
     bzero( &newtio, sizeof( newtio ) ); 

    
     newtio.c_cflag  |=  CLOCAL | CREAD;  
     newtio.c_cflag &= ~CSIZE;  

     switch( nBits ) 
     { 
     case 7: 
      newtio.c_cflag |= CS7; 
      break; 
     case 8: 
      newtio.c_cflag |= CS8; 
      break; 
     } 

     switch( nEvent ) 
     { 
     case 'o':
     case 'O': 
      newtio.c_cflag |= PARENB; 
      newtio.c_cflag |= PARODD; 
      newtio.c_iflag |= (INPCK | ISTRIP); 
      break; 
     case 'e':
     case 'E': 
      newtio.c_iflag |= (INPCK | ISTRIP); 
      newtio.c_cflag |= PARENB; 
      newtio.c_cflag &= ~PARODD; 
      break;
     case 'n':
     case 'N': 
      newtio.c_cflag &= ~PARENB; 
      break;
     default:
      break;
     } 

switch( nSpeed ) 
     { 
     case 2400: 
      cfsetispeed(&newtio, B2400); 
      cfsetospeed(&newtio, B2400); 
      break; 
     case 4800: 
      cfsetispeed(&newtio, B4800); 
      cfsetospeed(&newtio, B4800); 
      break; 
     case 9600: 
      cfsetispeed(&newtio, B9600); 
      cfsetospeed(&newtio, B9600); 
      break; 
     case 115200: 
      cfsetispeed(&newtio, B115200); 
      cfsetospeed(&newtio, B115200); 
      break; 
     case 460800: 
      cfsetispeed(&newtio, B460800); 
      cfsetospeed(&newtio, B460800); 
      break; 
     default: 
      cfsetispeed(&newtio, B9600); 
      cfsetospeed(&newtio, B9600); 
     break; 
     } 

     if( nStop == 1 ) 
      newtio.c_cflag &=  ~CSTOPB; 
     else if ( nStop == 2 ) 
      newtio.c_cflag |=  CSTOPB; 

     newtio.c_cc[VTIME]  = 0; 
     newtio.c_cc[VMIN] = 0; 

     tcflush(fd,TCIFLUSH); 

if((tcsetattr(fd,TCSANOW,&newtio))!=0) 
     { 
      perror("com set error"); 
      return -1; 
     } 
     printf("set done!\n"); 
     return 0; 
} 

double thetaVOffset=0;
double* Method1(double detaXa,double detaYa,double pos1a,double pos2a,double kVxa,double kVya )
{
	double inc1=0,inc2=0;
	double detaVx,detaVy;
	double thetaV;
	double posAbsDiff=fabs(pos1a-pos2a);
    int i=0;
    if(detaXa ==0&&detaYa == 0){
        i=10;
    }
             ddetaXa0 = ddetaXa1;
             ddetaXa1 = ddetaXa2;
             ddetaXa2 = ddetaXa3;
             ddetaXa3 = ddetaXa4;
             ddetaXa4 = ddetaXa5;
             ddetaXa5 = ddetaXa6;
             ddetaXa6 = ddetaXa7;
             ddetaXa7 = ddetaXa8;
             ddetaXa8 = ddetaXa9;
             ddetaXa9 = detaXa;
             ddetaYa0 = ddetaYa1;
             ddetaYa1 = ddetaYa2;
             ddetaYa2 = ddetaYa3;
             ddetaYa3 = ddetaYa4;
             ddetaYa4 = ddetaYa5;
             ddetaYa5 = ddetaYa6;
             ddetaYa6 = ddetaYa7;
             ddetaYa7 = ddetaYa8;
             ddetaYa8 = ddetaYa9;
             ddetaYa9 = detaYa;
            if(i>1){
                detaXa =KI0*detaXa;
                detaYa =KI0*detaYa;
                i--;
            }

            else{
                if(detaXa>0&&detaYa>0){
                    detaXa =kp*detaXa+ (ddetaXa0+ddetaXa1+ddetaXa2+ddetaXa3+ddetaXa4+ddetaXa5+ddetaXa6+ddetaXa7+ddetaXa8+ddetaXa9)*ki+kd*(detaXa-ddetaXa0);
                    detaYa = kp*detaYa+ (ddetaYa0+ddetaYa1+ddetaYa2+ddetaYa3+ddetaYa4+ddetaYa5+ddetaYa6+ddetaYa7+ddetaYa8+ddetaYa9)*ki+kd*(detaYa-ddetaYa0);
                }
                else if(detaXa<0&&detaYa<0){
                    detaXa =kp*detaXa+kd*(detaXa-ddetaXa0);
                    detaYa = kp*detaYa+kd*(detaYa-ddetaYa0);
                }
                else if(detaXa<0&&detaYa>0){
                    detaXa =kp*detaXa+kd*(detaXa-ddetaXa0);
                    detaYa = kp*detaYa+ (ddetaYa0+ddetaYa1+ddetaYa2+ddetaYa3+ddetaYa4+ddetaYa5+ddetaYa6+ddetaYa7+ddetaYa8+ddetaYa9)*ki+kd*(detaYa-ddetaYa0);
                }
                else{
                    detaXa =kp*detaXa+ (ddetaXa0+ddetaXa1+ddetaXa2+ddetaXa3+ddetaXa4+ddetaXa5+ddetaXa6+ddetaXa7+ddetaXa8+ddetaXa9)*ki+kd*(detaXa-ddetaXa0);
                    detaYa = kp*detaYa+kd*(detaYa-ddetaYa0);
                }
             }
            
	if( posAbsDiff==180){
		thetaV=-1;
		detaVx=0;
		detaVy=0;
		inc1=0;//00000000000000000000000000000000000000000
		inc2=0;
	}
	else{
		if (posAbsDiff > 180){thetaV = mod((pos1a + pos2a) / 2 + 180 + thetaVOffset, 360);}
		else{thetaV = mod((pos1a + pos2a) / 2 + thetaVOffset, 360);}
       
		double thetaCCD;
		if (detaXa == 0 && detaYa == 0)
		{thetaCCD = -1;}
		else if (detaXa == 0){thetaCCD = 0.5 * PI * (2 - fabs(detaYa) / detaYa);}
		else if (detaYa == 0){thetaCCD = 0.5 * PI * (1 - fabs(detaXa) / detaXa);}
        else if (detaXa > 0 && detaYa > 0){thetaCCD = atan(detaYa / detaXa);}
        else if (detaXa > 0 && detaYa < 0){thetaCCD = atan(detaYa / detaXa) + 2 * PI;}
        else{thetaCCD = atan(detaYa / detaXa) + PI;}
		double absPIxel = sqrt(pow(detaXa, 2) + pow(detaYa, 2));
		if (thetaCCD != -1){
			double thetaVx = mod(thetaV - 90, 360);
			thetaCCD = thetaCCD * 180 / PI;
			double thetaVCCD = mod(thetaCCD - thetaVx, 360); 
            thetaVCCD = thetaVCCD * PI / 180;
            detaVx = absPIxel * cos(thetaVCCD);
            detaVy = absPIxel * sin(thetaVCCD);
            detaVx = detaVx / 1280 * 3.1826 * PI / 180;
            detaVy = detaVy / 1024 * 2.5461 * PI / 180;
		}
		else{ detaVx = 0;detaVy = 0;}
        double L1=10*PI/180;
		double alpha;
		if(posAbsDiff < 180){alpha=posAbsDiff;}
        else{alpha=360-posAbsDiff;}
        alpha=alpha*PI/180;
        L0=2*L1*cos(alpha/2);
        double detaVx1 = (L0*cos((pos1a+pos2a)/2)+0.85)/0.9 * fabs(detaVx);//000000000000000000000000000000000000000000000000
        double detaVy1 = (L0*sin((pos1a+pos2a)/2)+0.85)/0.9 * fabs(detaVy);//000000000000000000000000000000000000000000000000
        // double detaVx1 = (tan(L0/(2*L1)/9*PI)-tan(L0/(2*L1)/9*PI-0.01))/tan(0.01)*fabs(detaVx);//000000000000000000000000000000000000000000000000
        // double detaVy1 = (tan(L0/(2*L1)/9*PI)-tan(L0/(2*L1)/9*PI-0.01))/tan(0.01)*fabs(detaVy);//000000000000000000000000000000000000000000000000
        // double detaVx1 = fabs(detaVx);//000000000000000000000000000000000000000000000000
        // double detaVy1 = fabs(detaVy);//000000000000000000000000000000000000000000000000
   
        double thetaX = atan(detaVx1 / (L0 + detaVy1));
		double thetaY;
		double aaaa=(detaVy +L0)/ (2 * L1 * cos(thetaX));
		if(aaaa>1){aaaa=1;}
		if(aaaa<-1){aaaa=-1;}
        thetaY =fabs(alpha / 2 - acos(aaaa));
        double detaVx2 = fabs(thetaX / PI * 180);
        double detaVy2 = fabs(thetaY / PI * 180);

		if (detaVx >= 0){
			if (posAbsDiff < 180){
				if (detaVy >= 0){
					if (pos2a > pos1a){
						inc1 = -kVxa * detaVx2 + kVya * detaVy2;
						inc2 = -kVxa * detaVx2 - kVya * detaVy2;
					}
					else if (pos2a == pos1a){
						inc1 = -kVxa * detaVx2;
						inc2 = -kVxa * detaVx2;
					}
					else{
						inc1 = -kVxa * detaVx2 - kVya * detaVy2;
						inc2 = -kVxa * detaVx2 + kVya * detaVy2;
					}
				}
				else{
					if (pos2a >= pos1a){
						inc1 = -kVxa * detaVx2 - kVya * detaVy2;
						inc2 = -kVxa * detaVx2 + kVya * detaVy2;
					}
					else{
						inc1 = -kVxa * detaVx2 + kVya * detaVy2;
						inc2 = -kVxa * detaVx2 - kVya * detaVy2;
					}
				}
			}
			else{
				if (detaVy >= 0){
					if (pos2a > pos1a){
						inc1 = -kVxa * detaVx2 - kVya * detaVy2;
						inc2 = -kVxa * detaVx2 + kVya * detaVy2;
					}
					else{
						inc1 = -kVxa * detaVx2 + kVya * detaVy2;
						inc2 = -kVxa * detaVx2 - kVya * detaVy2;
					}
				}
				else{
					if (pos2a >= pos1a){
						inc1 = -kVxa * detaVx2 + kVya * detaVy2;
						inc2 = -kVxa * detaVx2 - kVya * detaVy2;
					}
					else{
						inc1 = -kVxa * detaVx2 - kVya * detaVy2;
						inc2 = -kVxa * detaVx2 + kVya * detaVy2;
					}
				}
			}
		}
		else{
			if (posAbsDiff < 180){
				if (detaVy >= 0){
					if (pos2a > pos1a){
						inc1 = kVxa * detaVx2 + kVya * detaVy2;
                        inc2 = kVxa * detaVx2 - kVya * detaVy2;
					}
					else if (pos2a == pos1a){
						inc1 = kVxa * detaVx2;
						inc2 = kVxa * detaVx2;
					}
					else{
						inc1 = kVxa * detaVx2 - kVya * detaVy2;
						inc2 = kVxa * detaVx2 + kVya * detaVy2;
					}
				}
				else{
					if (pos2a >= pos1a){
						inc1 = kVxa * detaVx2 - kVya * detaVy2;
						inc2 = kVxa * detaVx2 + kVya * detaVy2;
					}
					else{
						inc1 = kVxa * detaVx2 + kVya * detaVy2;
						inc2 = kVxa * detaVx2 - kVya * detaVy2;
					}
				}
			}
			else{
				if (detaVy >= 0){
					if (pos2a > pos1a){
						inc1 = kVxa * detaVx2 - kVya * detaVy2;
						inc2 = kVxa * detaVx2 + kVya * detaVy2;
					}
					else{
						inc1 = kVxa * detaVx2 + kVya * detaVy2;
						inc2 = kVxa * detaVx2 - kVya * detaVy2;
					}
				}
				else{
					if (pos2a >= pos1a){
						inc1 = kVxa * detaVx2 + kVya * detaVy2;
						inc2 = kVxa * detaVx2 - kVya * detaVy2;
					}
					else{
						inc1 = kVxa * detaVx2 - kVya * detaVy2;
						inc2 = kVxa * detaVx2 + kVya * detaVy2;
					}
				}
			}
		}
	}

	static double posInc[2];
	posInc[0]=inc1;
	posInc[1]=inc2;
	return posInc;
}



// double* InverseParaxial(double Ep, double Ap)
// {
//     double theta11, theta12, theta21, theta22;
//     double zero;

//     if (Ep == 0)
//     {
//         zero = 0;
//         theta11 = mod(Ap, 2 * PI);
//         theta12 = mod(Ap + PI, 2 * PI);
//         theta21 = mod(Ap + PI, 2 * PI);
//         theta22 = mod(Ap, 2 * PI);
//     }
//     else
//     {
//         zero = 1;
//         double E1 = 10 * PI / 180;
//         double E2 = 10 * PI / 180;
//         double CC = (pow(Ep, 2) + pow(E1, 2) - pow(E2, 2)) / (2 * Ep * E1);
//         double CC1 = (pow(Ep, 2) + pow(E2, 2) - pow(E1, 2)) / (2 * Ep * E2);

//         double alpha, beta;
//         if (CC >= 1)
//         {
//             alpha = 0;
//         }
//         else
//         {
//             alpha = acos(CC);
//         }
//         if (CC1 >= 1)
//         {
//             beta = 0;
//         }
//         else
//         {
//             beta = acos(CC1);
//         }

//         theta11 = mod(Ap - alpha, 2 * PI);
//         theta12 = mod(Ap + beta, 2 * PI);
//         theta21 = mod(Ap + alpha, 2 * PI);
//         theta22 = mod(Ap - beta, 2 * PI);
//     }

//     theta11 = theta11 * 180 / PI;
//     theta12 = theta12 * 180 / PI;
//     theta21 = theta21 * 180 / PI;
//     theta22 = theta22 * 180 / PI;

//     double* inverseSolution;
//     inverseSolution[0] = theta11;
//     inverseSolution[1] = theta12;
//     inverseSolution[2] = theta21;
//     inverseSolution[3] = theta22;
//     inverseSolution[4] = zero;
//     return inverseSolution;
// }

// void * XIANCHENG6(void){
//     while(1){
//     }
// }



// void * XIANCHENG5(void){
//     char chaxunweizhi[4]={'p','f','b',0x0D};
//     while(1){
//         // write(fd4,chaxunweizhi,4);
//         // write(fd5,chaxunweizhi,4);
//         sleep(0.02);
//     }
// }


void * XIANCHENG4(void){

    while(trackStatus){
        if ( tuoBaLiangStatus)
                {
                        time_t stop4=clock();
                        write(fd4,chaxunweizhi,4);
                        write(fd5,chaxunweizhi,4);
                        double detaX1 = (double)detaX;
                        double detaY1 = (double)detaY;
                        tuoBaLiangStatus=0x00;
                        double pos1a = mod(pos1, 360);
                        double pos2a = mod(pos2, 360);
                        // printf("%f,%f\n" ,pos1a,pos2a);
                        if(detaX1==0){systimes=0;}
                        double x =PI/18 * (cos(pos1a*PI/180 )+cos(pos2a*PI/180 ))*1280 /3.1826 +detaX1;
                        double y =PI/18 * (sin(pos1a*PI/180 )+sin(pos2a*PI/180 ))* 1024 / 2.5461 +detaY1;
                        double *gj,*step;
                        step = Akima_forecast(systimes, x, systimes + 1, 0);
                        gj = Akima_forecast(systimes, y, systimes+1,1);
                        systimes ++ ;
                        // detaX1=step[1]-PI/18 * (cos(pos1a*PI/180 )+cos(pos2a*PI/180 ))*1280 /3.1826 ;
                        // detaY1=gj[1]-PI/18 * (sin(pos1a*PI/180 )+sin(pos2a*PI/180 ))* 1024 / 2.5461;
                        printf("x0=%f,y0=%f,x1=%f,y1=%f",x,y,step[1],gj[1]);
                        double * INC;
                        INC = Method1(detaX1, detaY1, pos1a, pos2a, kRotary, kAngle);
                        double cmdInc1 = INC[0];
                        double cmdInc2 = INC[1];
                        char* cm11;
                        cm11=ToString(cmdInc1);
                        char* cm1[strlen(cm11)];
                        strcpy(cm1,cm11);
                        char* cm21;
                        cm21=ToString(cmdInc2);
                        char* cm2[strlen(cm21)];
                        strcpy(cm2,cm21);
                        double v01 =fabs(cmdInc1*kvv+0.1);
                        char* v11 = ToString(v01);
                        char* v1[strlen(v11)];
                        strcpy(v1,v11);
                        double v02 =fabs(cmdInc2*kvv+0.1);
                        char* v22 = ToString(v02);
                        char* v2[strlen(v22)];
                        strcpy(v2,v22);
                        // printf("v1=%s,v2=%s\n",v1,v2);
                        char* mover ="moveinc ";
                        char* kongge = " ";
                        char* wu =" 1";
                        char* huiche[1]={0x0d};
                        char* pos1Cmd[strlen(mover)+strlen(cm1)+strlen(kongge)+strlen(v1)+strlen(wu)+strlen(huiche)+1];
                        char* pos2Cmd[strlen(mover)+strlen(cm2)+strlen(kongge)+strlen(v2)+strlen(wu)+strlen(huiche)+1];
                        strcpy(pos1Cmd,mover);
                        strcat(pos1Cmd,cm1);
                        strcat(pos1Cmd,kongge);
                        strcat(pos1Cmd,v1);
                        strcat(pos1Cmd,wu);
                        strcpy(pos2Cmd,mover);
                        strcat(pos2Cmd,cm2);
                        strcat(pos2Cmd,kongge);
                        strcat(pos2Cmd,v2);
                        strcat(pos2Cmd,wu);
                        strcat(pos1Cmd,huiche);
                        strcat(pos2Cmd,huiche);
                        // printf("pos1Cmd=%s\n",pos1Cmd);
                        // printf("pos2Cmd=%s\n",pos2Cmd);
                        write(fd4,pos1Cmd, strlen(pos1Cmd));
                        write(fd5,pos2Cmd, strlen(pos2Cmd));
                        stop3=clock();
                        time_doing_2 =(double)(stop3-stop4)/CLOCKS_PER_SEC;
                        // printf("using time2 = %f s\n",time_doing_2);
                    }
                else { continue; }
    }
    
}

void * XIANCHENG1(void){
    int counter=0;
    int con =0;
    for(;;){
        // usleep(1000);
        get6=read(fd6,buff6,1);
        if(get6>0){
            if(buff6[0]==0x43&&counter==0){
                TC[0]=0x43;
                counter=1;
            }
            else if(counter>0){
                TC[counter]=buff6[0];
                counter++;
            }
            }
             if(counter==13){
                unsigned char checksunm_0=ChecksumBit8(TC,12);
                if(TC[12]==checksunm_0&&TC[1]==0x63&&TC[2]==0x0c){
                switch (TC[3])
                {
                    case 0x01:
                        /* code */
                    break;
                    case 0x02:
                        if (TC[4]==0x01)
                       {
                            trackStatus=0x01;
                            int ret;
                            write(fd4,enable,3);
                            write(fd5,enable,3);
                            pthread_t Pthread[NUM_THREADS];
                            ret = pthread_create(&Pthread[3],NULL,(void *)XIANCHENG4,NULL);
                            stop1=clock();
                            time_doing_1 =(double)(stop1)/CLOCKS_PER_SEC;
                            // printf("using time0 = %f s\n",time_doing_1);
                            if(0!=ret){
                               printf("Error:defaut!");
                            }

                        }
                       else{trackStatus=0x00;write(fd4,unable,2);write(fd5,unable,2);}
                         break;
                    case 0x03:
                        /* code */
                        break;
                    case 0x04:
                        detaX=CharArrayToInt(TC[4],TC[5],TC[6],TC[7]);
                        detaY=CharArrayToInt(TC[8],TC[9],TC[10],TC[11]);
                        stop2 =clock();
                        time_doing_1 =(double)(stop2 )/CLOCKS_PER_SEC;
                        printf("using time1 = %f s\n",time_doing_1);
                        tuoBaLiangStatus=0x01;
                        break;
                    case 0x05:
                        /* code */
                        break;
                    default:
                        break;
                    }
                    con+=1;
                    counter=0;
                    // printf("%0x,%0x,%0x,%0x,%0x,%0x,%0x,%0x,%0x,%0x,%0x,%0x,%0x\n",TC[0],TC[1],TC[2],TC[3],TC[4],TC[5],TC[6],TC[7],TC[8],TC[9],TC[10],TC[11],TC[12]);
                    // printf("count=%d detaX=%d detaY=%d track=%d\n",con,detaX,detaY,trackStatus);
                    // printf("-----from xiancheng1\n");
                }
            }
        }
}

void * XIANCHENG2(void){
    int counter=0;
    int con =0;
    for(;;){
        // usleep(1000);
        get4=read(fd4,buff4,1);
        if(get4>0){
            if(buff4[0]=='p'&&counter==0){
                TC_1[0]='p';
                counter=1;
            }
            else if(counter>0){
                TC_1[counter]=buff4[0];
                counter++;
            }
            if(buff4[0]=='['){
                if(TC_1[1]=='f'&&TC_1[2]=='b'){
                    char tc1[counter-7];
                    for (int i = 0; i < counter-7; i++)
                    {
                        tc1[i]=TC_1[i+5];
                    }
                    pos1=atof(tc1);
                    counter=0;
                    con ++;
                    // printf("%s\n",TC_1);
                    // printf("pos1=%lf------con=%d-----from xiancheng2\n",pos1,con);
                }
            }
            }
        }

}

void * XIANCHENG3(void){
    int counter=0;
    int con =0;
    for(;;){
        // usleep(1000);
        get5=read(fd5,buff5,1);
        if(get5>0){
            if(buff5[0]=='p'&&counter==0){
                TC_2[0]='p';
                counter=1;
            }
            else if(counter>0){
                TC_2[counter]=buff5[0];
                counter++;
            }
            if(buff5[0]=='['){
                if(TC_2[1]=='f'&&TC_2[2]=='b'){
                    char tc2[counter-7];
                    for (int i = 0; i < counter-7; i++)
                    {
                        tc2[i]=TC_2[i+5];
                    }
                    pos2=atof(tc2);
                    counter=0;
                    con ++;
                    // printf("%s\n",TC_2);
                    // printf("pos2=%lf------con=%d-----from xiancheng2\n",pos2,con);
                }
            }
            }
        }

}



 //mian######################################################################

int main(void)
{
    start=clock();
    fd4 = open_port(fd4,4);
    fd5 = open_port(fd5,5);
    fd6 = open_port(fd6,6);
    set_opt(fd4,115200,8,'N',1);
    set_opt(fd5,115200,8,'N',1);
    set_opt(fd6,115200,8,'N',1);
    int  counter=0;
    int con=0;
    pthread_t Pthread[NUM_THREADS];
    int ret;
    ret = pthread_create(&Pthread[0],NULL,(void *)XIANCHENG1,NULL);
    ret = pthread_create(&Pthread[1],NULL,(void *)XIANCHENG2,NULL);
    // ret = pthread_create(&Pthread[4],NULL,(void *)XIANCHENG5,NULL);
    ret = pthread_create(&Pthread[2],NULL,(void *)XIANCHENG3,NULL);
    // ret = pthread_create(&Pthread[5],NULL,(void *)XIANCHENG6,NULL);


pthread_exit(NULL);
        close(fd4);
        close(fd5);
        close(fd6);
        return 0;
    }



   
