//#define __DOUBLE_PRECISION__
#include "../defs.h" // определения типов, инклюды

// константы расчета (чтобы долго не искать в main'е)
#define N 200
#define H_0 350
#define R 0.06647221276
#define alpha -0.5410520681
#define beta Pi/2
#define Eps 0.01
//#define Vo5_0 0.754
#define xi_inf 0.004973591970
#define P1 290.3761507
#define G 697.8648613
#define Ms 9.6
//#define u_0 0.134;
#define C 0.134;

int space_pressed = 0;
// глобальные переменные - координаты точек для рисования
int nez = 0;
int debug = 0, shot_finished;
quad h = 1.0/N;
int curve_point = 1, Nc, ic=0;
int cross_flag = NoCrossing, crtmp;
quad u_0=1,a=R,b=5*R, a_prev, b_prev, H_angle=beta;
quad rhoc = 0, tetac = 10, uc = u_0;
quad Vol1 = 0, Vol2 = 0, delta = 1;

FILE *f;
quad u[4*N];
quad rho[4*N];
quad teta[4*N];

//#include "../fs.h" // стандартные функции f1,f2,f3,cross
#include "../glfs.h" // функции рисования


quad H(quad r,quad z)
{
	quad tmp;
	quad A=-1.0*(R*R*R);
	tmp=9*A*A*z*z*r*r/(r*r+z*z)/(r*r+z*z)/(r*r+z*z)/(r*r+z*z)/(r*r+z*z)+
	(1+A*(1/__sqrtq((r*r+z*z)*(r*r+z*z)*(r*r+z*z))
	-3*z*z/__sqrtq((r*r+z*z)*(r*r+z*z)*(r*r+z*z)*(r*r+z*z)*(r*r+z*z))))*
	(1+A*(1/__sqrtq((r*r+z*z)*(r*r+z*z)*(r*r+z*z))
	-3*z*z/__sqrtq((r*r+z*z)*(r*r+z*z)*(r*r+z*z)*(r*r+z*z)*(r*r+z*z))));
	tmp=__sqrtq(tmp);
	//printf("H=%f\n",tmp);
	return tmp;
	//printf("H=%f, A=%f\n",tmp,A);
}

quad L(quad r,quad z,quad H0)
{
	quad xi = xi_inf*H0;
	quad H_nd = H(r,z);
	quad pp = (__expq(xi*H_nd)+__expq(-xi*H_nd))/(__expq(xi*H_nd)-__expq(-xi*H_nd))-1/(xi*H_nd);
	//printf("P=%f\n",pp);
	return Ms*pp;
}

quad P(quad r,quad z,quad H0)
{
	quad xi = xi_inf*H0;
	quad H_nd = H(r,z);
	quad pp = __logq(__expq(xi*H_nd)/2-__expq(-xi*H_nd)/2)-__logq(xi*H_nd);
	//printf("P=%f\n",pp);
	return pp;
}

quad f3(quad r, quad u,quad C, quad H0, quad theta) {
	quad ff, k;
	k=-G*u+P1*P(r,u,H0);
	//printf("-Ah+P=%f,C-Ah+P=%f  ",k,k+C);
	if(r==0) {
		ff=(C+k)*0.5;
		return ff;
	}
	ff=(C+k)-1.0*__sinq(theta)/r;
	return ff;
}

quad f1(quad y) {
	return __cosq(y);
}
quad f2(quad y) {
	return -__sinq(y);
}

// CHECKING: CURVE CROSS HIMSELF AND DIRECTION OF THE TURN //
int cross(quad *x,quad *y,int n)
{int j;
quad f1,f2,g1,g2,k;
for(j=2;j<n-1;j++)
  {
  f1=(x[j]-x[n-1])/(x[n]-x[n-1])-(y[j]-y[n-1])/(y[n]-y[n-1]);
  f2=(x[j-1]-x[n-1])/(x[n]-x[n-1])-(y[j-1]-y[n-1])/(y[n]-y[n-1]);
  g1=(x[n]-x[j-1])/(x[j]-x[j-1])-(y[n]-y[j-1])/(y[j]-y[j-1]);
  g2=(x[n-1]-x[j-1])/(x[j]-x[j-1])-(y[n-1]-y[j-1])/(y[j]-y[j-1]);
  k=(x[j]-x[n])*(y[n-1]-y[j])-(y[j]-y[n])*(x[n-1]-x[j]);
  if((f1*f2<0)&&(g2*g1<0))
	 {if(k>0) return CounterClockWise;
	 if(k<0) return ClockWise;
	 }
  }
return NoCrossing;
}

char calc_step(int i, int flag, quad C){ 
/*	
	returns new flag using old cross_flag
	0: no crossing
	1, -1: self-crossing clockwise or counterclockwise
	10: crossover of the curve and sphere
	-10: curve is going away
*/
  	quad k11, k12, k13, k21, k22, k23;

  	k11=h*f1(teta[i-1]);
  	k12=h*f2(teta[i-1]);
  	k13=h*f3(rho[i-1], u[i-1], C, H_0, teta[i-1]);
  	k21=h*f1(teta[i-1]+k13);
  	k22=h*f2(teta[i-1]+k13);
  	k23=h*f3(rho[i-1]+k11,u[i-1]+k12,C,H_0,teta[i-1]+k13);
  	rho[i]=rho[i-1]+(k11+k21)*0.5;
  	u[i]=u[i-1]+(k12+k22)*0.5;
  	teta[i]=teta[i-1]+(k13+k23)*0.5;
  	
  	// CHECKING CROSS //
   	if((rho[i]<0)&&(u[i]<u[0]))
   		flag = ClockWise;
   	if((rho[i]<0)&&(u[i]>u[0]))
   		flag = CounterClockWise;
	if(i>2){
		crtmp = cross(rho,u,i);
		if(
			((crtmp == CounterClockWise) && (flag == ClockWise))
			||
			((crtmp == ClockWise) && (flag == CounterClockWise))
		)
		flag = NoCrossing;
		else if(crtmp != NoCrossing)
			flag = crtmp;
	}
  	
  	if(rho[i]*rho[i]+u[i]*u[i]<R*R){
  		flag = SphereCrossing;
  		//printf("curve&sphere: rho=%f, u=%f\n",(double)rho[i], (double)u[i]);
  	}	
  	if(i>3*N)
  		flag = GoesToInfty;

return flag;
} // end of calc_step  

void one_point_calc(){
	cross_flag = calc_step(curve_point, cross_flag, C);
	tetac = teta[curve_point];
	rhoc = rho[curve_point];
	uc = u[curve_point];
	curve_point++;
}

quad shoot(quad u_0_old){
	printf("u_0=%f,teta=%f, rhoc=%f\n",(double)u_0,(double)(tetac+Pi/2),(double)rhoc);
	if(nez>0){
		a = a_prev; b = b_prev;
		u_0 = (a+b)/2;
		if(nez<5)
			b=u_0;
		else 
			a=u_0;
		u_0 = (a+b)/2;
	}

	a_prev = a; b_prev = b;

	while((rhoc<1)&&(cross_flag == NoCrossing))
		one_point_calc();

	if(nez == 0) {
		if(cross_flag == SphereCrossing){
			printf("Drop\n");
			b=u_0_old;
		}
		if(cross_flag == ClockWise)
			b=u_0_old;
		if(cross_flag == CounterClockWise)
			a=u_0_old;
		if(cross_flag == NoCrossing){
			if(tetac-alpha<-Eps)
				a=u_0_old;
			if(tetac-alpha>Eps)
				b=u_0_old;
		}
		if((cross_flag == GoesToInfty)&&(uc < u[0]))
			a=u_0_old;
		if((cross_flag == GoesToInfty)&&(uc > u[0]))
			b=u_0_old;
	}
	cross_flag = NoCrossing;
	printf("new u_0=%f\n",(double)(a+b)/2);
	printf("delta_u_0=%e\n",(double)(u_0_old-(a+b)/2));
	return (a+b)/2;
}

double MF(double x, double y, double H0){
	double rhof, P0, ff;
	quad qx, qy;
	qx = (quad)x; qy = (quad)y;
	ff=0;
	
	P0=2974.325581;
	ff=H0*H0*pow((double)H(qx,qy),2)/(8*Pi)-P0*(double)P(qx,qy,(quad)H0)+(double)L(qx,qy,(quad)H0)*H0*(double)H(qx,qy)/4/Pi;
	return ff;
}



void change_debug_mode(){
	if(debug)
		debug = 0;
	else
		debug = 1;
}


void quit_and_save(){
	int i;
	quad Vol=0;
	printf("Quitting: TETAsm=%f\n",(double)tetac);
	printf("\nkey pressed:%d times\n\n", space_pressed);
	printf("h0=%f, C=%f\n",(double)u[0],(double)C);

	
	// CALCULATIG OF THE VOLUME //
	for(i=1;i<Nc;i++){
		Vol=Vol+0.5*Pi*(u[i-1]+R+u[i]+R)*(rho[i]*rho[i]-rho[i-1]*rho[i-1]);
	}
	Vol=Vol-4*Pi*R*R*R/3;
	for(i=0;i<Nc+1;i++)
{ if(rho[i]*rho[i]+u[i]*u[i]-R*R<delta)
  {delta=rho[i]*rho[i]+u[i]*u[i]-R*R;ic=i;}
}
printf("free_mag_energy=%f; gravity energy=%f; surftension energy=%f\n",fme, fge, fse);
	printf("V0l=%f, rho0=%f, Volleft=%f, Volright=%f \n",(double)Vol,(double)rho[ic],(double)Vol1,(double)Vol2);
	fprintf(f,"C=%f; V=%f; h0=%f; H0=%f;\n",(double)C,(double)Vol,(double)u[0],(double)H_0);
	//fprintf(f,"free_mag_energy=%f; gravity energy=%f; surftension energy=%f\n",fme, fge, fse);
	for(i=0;i<Nc;i++)
	fprintf(f,"[%f,%f], ",(double)rho[i],(double)u[i]);
	fclose(f);
	scanf("%d",&i);
	exit(1);
}

void keyWindow (unsigned char k, int x, int y) {
	
	quad u_0_old = u_0;
	space_pressed++;
	nez = 0;
	switch (k) {
    	case 0x1B : exit(0);
	    case '1': nez = 1; break;
    	case '2': nez = 10; break;
		case '3': nez = 666; break;
	    case '9': change_debug_mode(); break;
	}
	if(
		((fabs((double)(tetac-alpha))<Eps))
		||
		(nez == 666)
		)
		quit_and_save();
  	if(debug==1){
  		one_point_calc();
  		printf("rho=%f, flag=%d\n",(double)(rhoc), cross_flag);
  		if((rhoc>=1)||(cross_flag!=NoCrossing)){
  			u_0 = shoot(u_0);
			if(u_0-u_0_old==0)
				printf("more precision needed\n");
			shot_finished = 1;
  		}
	} else {
		rhoc=rho[0];
		uc=u_0;
		u_0 = shoot(u_0);
		if(u_0-u_0_old==0)
			printf("more precision needed\na-C=%e, b-c=%e\n",(double)(a-C),(double)(b-C));
		shot_finished = 1;
	}
 	drawWindow();
	if(shot_finished){
		Nc = curve_point;
 		curve_point	= 1;
	}
}

// инициализация графики и начальные данные 
void userInterfaceInit(int argc, char **argv){
	// инициализация графики
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	
	// создание окна
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(900, 400);
	glutCreateWindow("Sphere");

	// настройка обработчиков событий
	glutDisplayFunc(drawWindow);
	glutReshapeFunc(reshapeWindow);
	glutKeyboardFunc(keyWindow);

	// цикл ожидания нажатия клавиши (keyWindow)
	glutMainLoop();
	
}

            
int main(int argc, char **argv)
{
	f=fopen("sphere.txt","w");
	//printf("xi_0=%f, H_0=%f, Tetasm=%f\n",xi_inf,H_0,-Pi/2);
	
	u[0]=u_0;
	rho[0] = 0;
	tetac = 1000;
	teta[0] = 0;

	userInterfaceInit(argc, argv);
	
	return 0;
}
