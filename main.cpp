#include <stdio.h>
#include <math.h>
#include <malloc.h>
#define gamma 1.4

#define NoError 0x00
#define ErrorMemory 0x01


typedef double hlleReal;
typedef int hlleInt;
typedef unsigned int hlleUInt;

typedef struct
{
    hlleReal F1;//rho*u
    hlleReal F2;//rhi*u^2+p
    hlleReal F3;//(e+p)*u
}hlleFlux;

typedef struct
{
    hlleReal U1;//rho
    hlleReal U2;//rho*u
    hlleReal U3;//e=rho*eps+rho*u^2/2 eps=gamma/(gamma-1)p/rho
    hlleReal u;
    hlleReal p;
    hlleReal eps;
}hlleConservative;

FILE *outRho,*outP,*outU,*outRhoU,*outE,*outEps;
hlleUInt T;//время
hlleUInt N;//ось Х

hlleFlux F(hlleConservative U);
hlleReal a(const hlleReal x,const hlleReal y);
hlleConservative hlleUpdateConservative(hlleConservative U);
hlleFlux hlleCalculationFlux(const hlleConservative U_Left,const hlleConservative U_Right);
void hlleInitialConditions(hlleConservative *U_L,hlleConservative *U_R);
void hlleOpenFiles(void);
void hlleCloseFiles(void);
void hlleWriteFiles(hlleConservative U);
void hlleWriteFiles(void);


int main(void)
{
    hlleFlux Flux,FluxPrevious;

    hlleConservative U_L,U_R;
    hlleConservative *U;

    hlleOpenFiles();
    hlleInitialConditions(&U_L,&U_R);
    U=(hlleConservative*) malloc(N*sizeof(hlleConservative));

    if(U==NULL)//
    {
        printf("The large size of the array");
        getchar();
        return ErrorMemory;
    }

    //U(x,0)=U_L(x) , x<X_0 or U_R(x), x>X_0
    for(hlleUInt i=0;i<N;++i)
    {

        if((double)i/N<0.5)
        {
            U[i]=U_L;
        }
        else
        {
            U[i]=U_R;

        }
    }

    for(hlleUInt n=0;n<T;++n)
    {
        FluxPrevious.F1=U[0].U1*U[0].u;
        FluxPrevious.F2=U[0].U1*U[0].u*U[0].u+U[0].p;
        FluxPrevious.F3=(U[0].U3+U[0].p)*U[0].u;
        //Flux=hlleCalculationFlux(U[n][0],U[n][1]);
        for(hlleUInt i=0;i<N-1;++i)
        {
            hlleWriteFiles(U[i]);//записаь значений в t=0;
            Flux=hlleCalculationFlux(U[i],U[i+1]);
            U[i].U1=U[i].U1-(N*1.0)/T*(Flux.F1-FluxPrevious.F1);
            U[i].U2=U[i].U2-(N*1.0)/T*(Flux.F2-FluxPrevious.F2);
            U[i].U3=U[i].U3-(N*1.0)/T*(Flux.F3-FluxPrevious.F3);
            U[i]=hlleUpdateConservative(U[i]);
            FluxPrevious=Flux;
        }
        U[N-1]=U_R;
        hlleWriteFiles();

    }
    hlleCloseFiles();
    free(U);
    return NoError;
}

hlleFlux F(hlleConservative U)
{
    hlleFlux F;
    F.F1=U.U1*U.u;//rho*u
    F.F2=U.U1*U.u*U.u+U.p;//rho*u^2+p
    F.F3=(U.U3+U.p)*U.u;//(e+p)*u
    return  F;
}

hlleReal a(const hlleReal x,const hlleReal y)// a(x,y)--max, -a(-x,-y)--min;
{
    if((x<0.0)&&(y<0.0))
    {
        return 0.0;
    }
    else
    {
        if(x>y)
        {
            return x;
        }
        else
        {
            return y;
        }
    }
}


hlleConservative hlleUpdateConservative(hlleConservative U)
{
    //предиктор?(корректор)
    U.u=U.U2/U.U1;
    U.eps=U.U3/U.U1-U.u*U.u/2;
    U.p=U.eps*U.U1*(gamma-1)/gamma;
    return U;
}

hlleFlux hlleCalculationFlux(const hlleConservative U_Left,const hlleConservative U_Right)
{
    hlleFlux returnFlux;
    hlleReal cRight,cLeft;
    hlleReal aPlus,aMinus;
    hlleReal uTilda,rhoTilda,cTilda,epsTilda;

    //Roe-average урседенние по Роу
    rhoTilda=sqrt(U_Left.U1*U_Right.U1);

    uTilda=sqrt(U_Left.U1)*U_Left.u+sqrt(U_Right.U1)*U_Right.u;
    uTilda=uTilda/(sqrt(U_Left.U1)+sqrt(U_Right.U1));

    epsTilda=sqrt(U_Left.U1)*U_Left.eps+sqrt(U_Right.U1)*U_Right.eps;
    epsTilda=epsTilda/(sqrt(U_Left.U1)+sqrt(U_Right.U1));

    cTilda=sqrt((gamma-1)*epsTilda)/rhoTilda;
    //

    //Скорость звука
//    cRight=sqrt(gamma*U_Left.p/U_Left.U1);
//    cLeft=sqrt(gamma*U_Right.p/U_Right.U1);
    cRight=sqrt((gamma-1)*U_Right.eps)/U_Right.U1;
    cLeft=sqrt((gamma-1)*U_Left.eps)/U_Left.U1;

    aPlus=a(U_Right.u+cRight,uTilda+cTilda);
    aMinus=-a(-(U_Left.u-cLeft),-(uTilda-cTilda));
    returnFlux.F1=aPlus*F(U_Left).F1-aMinus*F(U_Right).F1+aPlus*aMinus*(U_Right.U1-U_Left.U1);
    returnFlux.F2=aPlus*F(U_Left).F2-aMinus*F(U_Right).F2+aPlus*aMinus*(U_Right.U2-U_Left.U2);
    returnFlux.F3=aPlus*F(U_Left).F3-aMinus*F(U_Right).F3+aPlus*aMinus*(U_Right.U3-U_Left.U3);

    returnFlux.F1=returnFlux.F1/(aPlus-aMinus);
    returnFlux.F2=returnFlux.F2/(aPlus-aMinus);
    returnFlux.F3=returnFlux.F3/(aPlus-aMinus);

    return returnFlux;
}
void hlleInitialConditions(hlleConservative *U_L,hlleConservative *U_R)
{
    FILE *IN;
    IN=fopen("InitialConditions.dat","r");
    //Read N and T
    if(fscanf(IN,"N=%d;\n",&N))
    {
         printf("N=%d;\n",N);
    }
    else
    {
        printf("data from the file have not been read");
    }
    if(fscanf(IN,"T=%d;\n\n",&T))
    {
         printf("T=%d;\n\n",T);
    }
    else
    {
        printf("data from the file have not been read");
    }

    //Read left conditions
    if(fscanf(IN,"rhoLeft=%lf;\n",&U_L->U1))
    {
        printf("rhoLeft=%lf;\n",U_L->U1);
    }
    else
    {
        printf("data from the file have not been read");
    }
    if(fscanf(IN,"velocityLeft=%lf;\n",&U_L->u))
    {
        printf("velocityLeft=%lf;\n",U_L->u);
    }
    else
    {
        printf("data from the file have not been read");
    }
    if(fscanf(IN,"pressureLeft=%lf;\n",&U_L->p))
    {

        printf("pressureLeft=%lf;\n\n",U_L->p);
        U_L->U2=U_L->U1*U_L->u;
        U_L->eps=gamma/(gamma-1)*U_L->p/U_L->U1;
        U_L->U3=U_L->U1*U_L->eps+U_L->U1*U_L->u*U_L->u/2;
    }
    else
    {
        printf("data from the file have not been read");
    }

    //Read right conditions
    if(fscanf(IN,"rhoRight=%lf;\n",&U_R->U1))
    {
        printf("rhoLeft=%lf;\n",U_R->U1);
    }
    else
    {
        printf("data from the file have not been read");
    }
    if(fscanf(IN,"velocityRight=%lf;\n",&U_R->u))
    {
        printf("velocityLeft=%lf;\n",U_R->u);
    }
    else
    {
        printf("data from the file have not been read");
    }
    if(fscanf(IN,"pressureRight=%lf;\n",&U_R->p))
    {

        printf("pressureLeft=%lf;\n\n",U_R->p);
        U_R->U2=U_R->U1*U_R->u;
        U_R->eps=gamma/(gamma-1)*U_R->p/U_R->U1;
        U_R->U3=U_R->U1*U_R->eps+U_R->U1*U_R->u*U_R->u/2;
    }
    else
    {
        printf("data from the file have not been read");
    }
    fclose(IN);
}
void hlleOpenFiles(void)
{
    outRho=fopen("rho.dat","w+");
    outP=fopen("p.dat","w+");
    outU=fopen("u.dat","w+");
    outRhoU=fopen("rhoU.dat","w+");
    outE=fopen("E.dat","w+");
    outEps=fopen("Eps.dat","w+");
}
void hlleCloseFiles(void)
{
    fclose(outRho);
    fclose(outP);
    fclose(outU);
    fclose(outRhoU);
    fclose(outE);
    fclose(outEps);
}
void hlleWriteFiles(hlleConservative U)
{
    fprintf(outRho,"%lf\t",U.U1);
    fprintf(outP,"%lf\t",U.p);
    fprintf(outU,"%lf\t",U.u);
    fprintf(outRhoU,"%lf\t",U.U2);
    fprintf(outE,"%lf\t",U.U3);
    fprintf(outEps,"%lf\t",U.eps);
}
void hlleWriteFiles(void)
{
    fprintf(outRho,"\n");
    fprintf(outP,"\n");
    fprintf(outU,"\n");
    fprintf(outRhoU,"\n");
    fprintf(outE,"\n");
    fprintf(outEps,"\n");
}
