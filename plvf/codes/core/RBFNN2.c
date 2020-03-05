         /* p: the group number of train data
   m: the number of input layer
   n: the number of middle layer
   q: the number of output layer  */                    
#include "stdlib.h"
#include "stdio.h"
#include "conio.h"
#include "math.h"

//从文件中读取数据：1 训练数据 2 权值数据  3测试数据  4中心数据 
void  readfile(x,w,newmid,p,m,q,n)
float *x,*w,*newmid;
int m,p,q,n;

{
FILE *fp;
float mid;
int i,nn;
printf("\n(1):Train Date File:traindata.dat\n(2):Weight Date File:w.dat\n(3):Test Data File:testdata.dat\n(4):Middle Point Data File:c.dat\n");
scanf("%d",&nn);

if(nn==1)   //训练数据
  {
  if ((fp=fopen("F:\\vcpp\\Test\\myinfuc\\traindata.dat","r"))==NULL)
      printf("\n The file cannot be opened!\n");
  else
     for(i=0;i<p*(m+q);i++)
     {
	  fscanf(fp,"%f",&mid);//从fp所指文件中读入一个浮点数放入变量mid中
	  x[i]=mid;
     }

     fclose(fp);
  }

if(nn==2)//权值数据
 {

 if ((fp=fopen("F:\\vcpp\\Test\\myinfuc\\w.dat","r"))==NULL)
      printf("cannot open");
 else
      for(i=0;i<(n+1)*q;i++)
       {
	   fscanf(fp,"%f",&mid);
	    w[i]=mid;
       }

       fclose(fp);
  }

if(nn==3)//测试数据
  {
   if((fp=fopen("testdata.dat","r"))==NULL)
       printf("cannot open");
   else
      for(i=0;i<p*(m+q);i++)
      {
	fscanf(fp,"%f",&mid);
	x[i]=mid;
      }

      fclose(fp);
  }

if(nn==4)//中心数据
  {
     if ((fp=fopen("c.dat","r"))==NULL)
      printf("cannot open");
      else
       for(i=0;i<p*n;i++)
      {
	 fscanf(fp,"%f",&mid);
	 newmid[i]=mid;
      }

      fclose(fp);
  }

}

//用k-means方法确定隐层节点
int lookfor(float *x,float *newmid,int *numberofmid,int p,int m,int q,int n)

{
int i,j,k,k2,i2;
int flag,f,cx=0,ss;
float *x1,*old;
float *sum,u;
double *l;
x1=(float*) calloc (p*(m+q+1),sizeof(float));
old=(float*) calloc (n*(m+q+1),sizeof(float));
l=(double*) calloc (n,sizeof(double));
sum=(float*) calloc (m+q,sizeof(float));

for(i=0;i<p;i++)
  {
      x1[i*(m+q+1)]=0;
      for(j=0;j<m+q;j++)
		x1[i*(m+q+1)+j+1]=x[i*(m+q)+j];
  }

for(i=0;i<n;i++)
  {
     old[i*(m+q+1)]=(float)i;
     for(j=1;j<m+q+1;j++)
	   old[i*(m+q+1)+j]=x1[i*(m+q+1)+j];
  }

do
{
 flag=1;
 for(i=0;i<p;i++)
   {
     f=0;
     for(j=0;j<n;j++)
      {
	l[j]=0;
	for(k=1;k<m+q+1;k++)
	    l[j]+=(x1[i*(m+q+1)+k]-old[j*(m+q+1)+k])*(x1[i*(m+q+1)+k]-old[j*(m+q+1)+k]);
	l[j]=sqrt(l[j]);
	if(l[f]>l[j])  f=j;
      }
     x1[i*(m+q+1)]=(float)f;
   }
 i2=0;
 for(i=0;i<n;i++)
  {
   newmid[i*(m+q+1)]=(float)i;
   
   for(j=0;j<m+q;j++)
	   sum[j]=0;
   ss=0;
   for (k=0;k<p;k++)
   {
	   if(x1[k*(m+q+1)]==i2)
	   {
		   ss++;
		   for (j=1;j<m+q+1;j++)
		   {
			   sum[j-1]+=x1[k*(m+q+1)+j];
		   }
	   }
   }
	   
  numberofmid[i2]=ss;
  if(ss!=0)
	for(j=1;j<m+q+1;j++)
	{
		sum[j-1]/=ss;
		newmid[i*(m+q+1)+j]=sum[j-1];
	}
  else
		{i=i-1;n-=1;}
  i2++;
 }
 n=i2;
  for(j=0;j<n;j++)
  {
	  if(numberofmid[j]==0)
	  {
		  for(k=j;k<n-1;k++)
		  {
			  old[(k+1)*(m+q+1)]=(float)(k-1);
			  for(k2=1;k2<m+q+1;k2++)
				  old[k*(m+q+1)+k2]=old[(k+1)*(m+q+1)+k2];
		  }
	  }
  }
  n=i;
  
  for(i=0;i<n;i++)

       for(j=1;j<m+q+1;j++)
	  {
	  u=(float)fabs((double)(newmid[i*(m+q+1)+j]-old[i*(m+q+1)+j]));
	  if(u>1e-6) flag=0;
	  old[i*(m+q+1)+j]=newmid[i*(m+q+1)+j];
	  }
   cx++;
 }
while(flag==0);
printf("\n");
printf("cx=%d",cx);
return(n);
}

//数据预处理函数
void process(x,p,m,q)
float *x;
int m,p,q;
{
int i,j;
float max,min;
for(i=0;i<m+q;i++)
  {
   max=x[i];
   min=x[i];
   for(j=0;j<p;j++)
      {
       if(x[j*(m+q)+i]>max)   max=x[j*(m+q)+i];
       if(x[j*(m+q)+i]<min)   min=x[j*(m+q)+i];
      }
   for(j=0;j<p;j++)
      x[j*(m+q)+i]=(x[j*(m+q)+i]-min)/(max-min);
  }
return;
}

//隐层输出计算函数
void midfun(x,newmid,t,p,m,q,n)
float *x,*newmid;
double *t;
int p,m,q,n;
{
int i,j,k;
for(i=0;i<p;i++)
   for(j=0;j<n;j++)
    {
    t[i*n+j]=0;
    for(k=0;k<m;k++)
	t[i*n+j]+=(double)((x[i*(m+q)+k]-newmid[j*(m+q+1)+k+1])*(x[i*(m+q)+k]-newmid[j*(m+q+1)+k+1]));
    t[i*n+j]=exp(-t[i*n+j]/2);
    }
return;
}


//初始化输出层权阈值
void initial(P,W,n,q)
float *P,*W;
int n,q;
{
float a,c;
int i,j;
printf("please input a:");
scanf("%f",&a);
printf("\n");
for(i=0;i<n+1;i++)
   for(j=0;j<n+1;j++)
       P[i*(n+1)+i]=a*a;

printf("please input c:");
scanf("%f",&c);
printf("\n");

for(i=0;i<n+1;i++)
   for(j=0;j<q;j++)
      W[i*q+j]=c;
}



//矩阵相乘函数
void matrixM(float *A,float *B,float *C,int n1,int n2,int n3)
//A: n1*n2
//B: n2*n3
//C: n1*n3
{
	int i,j,k;
	for(i=0;i<n1;i++)
		for(j=0;j<n3;j++)
			for(k=0;k<n2;k++)
				C[i*n3+j]+=A[i*n2+k]*B[k*n3+j];

}

//用递推最小二乘法来训练输出权值
void get_w(W,P,H,d,A,n,q)

float *W,*P,*H,d,A;
int n,q;
{
 int i,j;
 float c=0,*PH,*K,*HPH,*F,*HH,*HW,*KHW,*P2;
 KHW=(float*) calloc (n+1,sizeof(float));
 HW=(float*) calloc (1,sizeof(double));
 HH=(float*) calloc ((n+1)*(n+1),sizeof(float));
 F=(float*) calloc (n+1,sizeof(float));
 K=(float*) calloc (n+1,sizeof(float));
 HPH=(float*) calloc (1,sizeof(float));
 PH=(float*) calloc (n+1,sizeof(float));
 P2=(float*) calloc ((n+1)*(n+1),sizeof(float));
 
 matrixM(P,H,PH,n+1,n+1,q);
 matrixM(H,PH,HPH,q,n+1,q);
 *HPH=(float)(1/(*HPH+1/A));
 matrixM(PH,HPH,K,n+1,q,q);
 matrixM(H,W,HW,q,n+1,q);
 *HW=(float)(d-*HW);
 matrixM(K,HW,KHW,n+1,q,q);
 for(i=0;i<n+1;i++)
	W[i]+=KHW[i];
 matrixM(K,H,HH,n+1,q,n+1);
 for (i=0;i<n+1;i++)
 {
	 HH[i*(n+1)+i]-=1;
	 for (j=0;j<n+1;j++)
	 		 HH[i*(n+1)+j]*=(-1);
 }
 matrixM(HH,P,P2,n+1,n+1,n+1);
 for(i=0;i<n+1;i++)
	 for (j=0;j<n+1;j++)
		 P[i*(n+1)+j]=P2[i*(n+1)+j];
}
 
//测试网络函数
void test(W,t1,y,n,p)
float *W,*y;
double *t1;
int n,p;
{
int i,j;
for(i=0;i<p;i++)
   {
   y[i]=0;
   for(j=0;j<n+1;j++)
      y[i]+=W[j]*(float)(t1[i*(n+1)+j]);
   }
}


//主函数
#define  PI 3.141596

main()
{
float *x,*newmid,*y,*dd;
double  *t,*t1;
float *P,*W,*H;
int p=4,m=2,q=1,n=3,a,b;

float A,d,*e;
int i,j,k,*numberofmid;
long int l;
P=(float*) calloc ((n+1)*(n+1),sizeof(float));
W=(float*) calloc ((n+1)*q,sizeof(float));
x=(float*) calloc (p*(m+q),sizeof(float));
newmid=(float*) calloc (n*(m+1),sizeof(float));
t=(double*) calloc (p*n,sizeof(double));
t1=(double*) calloc (p*(n+1),sizeof(double));
H=(float*) calloc ((n+1)*q,sizeof(float));
y=(float*) calloc (p*q,sizeof(float));
e=(float*) calloc (p*q,sizeof(float));
dd=(float*) calloc (p*q,sizeof(float));
numberofmid=(int*) calloc (n,sizeof(int));
printf("Will you train now?\n 1.yes.     2.no.\n");
scanf("%d",&b);
printf("\n");

x[0]=0;x[1]=0;x[2]=0;
x[3]=0;x[4]=1;x[5]=1;
x[6]=1;x[7]=0;x[8]=1;
x[9]=1;x[10]=1;x[11]=0;

/*
for(i=0;i<p;i++)
	{
		x[i*(m+q)]=i*0.2f-1;
        x[i*(m+q)+1]=(float)(sin(2*PI*x[i*(m+q)])*exp(-x[i*(m+q)]));
	}
*/
//readfile(x,W,newmid,p,m,q,n);
for(i=0;i<p;i++)
   {
	for(j=0;j<m+q;j++)
		printf("x[%d][%d]=%1.1f ",i,j,x[i*(m+q)+j]);
	printf("\n"); 
   }
while(b==1)
{
b=0;
/*readfile(x,W,newmid,p,m,q,n);*/
/*process(x,p,m,q);*/
n= lookfor(x,newmid,numberofmid,p,m,q,n);
printf("网络隐层节点即中心点的坐标:\n");
for(i=0;i<n;i++)
{
	for(j=1;j<m+q;j++)
		printf("%1.1f  ",newmid[i*(m+q+1)+j]);
	printf("\n");

}
//writew(newmid,p,n);
midfun(x,newmid,t,p,m,q,n);
printf("please  input A: (Beteen 0 and 1) ");
scanf("%f",&A);
printf("\n");
initial(P,W,n,q);

for(i=0;i<p;i++)
  {
	for(j=0;j<n;j++)
		   t1[i*(n+1)+j]=t[i*n+j];
	t1[i*(n+1)+n]=1.0;
  }
for(i=0;i<p;i++)
	 dd[i]=x[i*(m+q)+m];

do
{
for(l=0;l<500;l++)
   for(i=0;i<p;i++)
     {
       for(j=0;j<n+1;j++)
	   H[j]=(float)(t1[i*(n+1)+j]);
	d=dd[i];
	get_w(W,P,H,d,A,n,q);
     }

for(i=0;i<p;i++)
{
  e[i]=0;
  for(k=0;k<p;k++)
     {
      d=dd[k];
      y[k]=0;
      for(j=0;j<n+1;j++)
	 {
	  H[j]=(float)(t1[k*(n+1)+j]);
	  y[k]+=W[j]*H[j];
	 }
      e[i]+=(y[k]-dd[k])*(y[k]-dd[k]);
     }
}
printf("\n\n\ne=%f\n",e[0]);
//getch();
printf("\n\n\nGo on training?\n 1.yes.  2.no.\n");
scanf("%d",&a);
printf("\n");

}while(a==1);

printf("输出层权值w[n+1][q]:\n");

for(i=0;i<n+1;i++)
    for(j=0;j<q;j++)
       printf("%f  ",W[i*q+j]);

printf("\n");

test(W,t1,y,n,p);
printf("网络期望输出d和测试输出y:\n");

   for(i=0;i<p;i++)
   {
	printf("d%d=%f   y%d=%f  ",i,x[i*(m+q)+m],i,y[i]);
	 printf("\n"); 
   }

	printf("\n");


}
}