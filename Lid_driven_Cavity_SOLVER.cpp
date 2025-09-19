#include<iostream>
#include<iomanip>
#include<cmath>
#define N 128
#define epsilon 1e-05
using namespace std;

void intial_conditions_for_cell_center_velocities_and_Pressure(double u[N+2][N+2], double v[N+2][N+2], double P[N+2][N+2])
{
    for(int i = 0; i < N+2; i++)
    {
        for(int j = 0; j < N+2; j++)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            P[i][j] = 0.0;
        }
    }
    for(int i = 1; i<=N; i++)
    {
        // CORRECTED: Changed initial lid velocity from 2.0 to 1.0
        u[N+1][i] = 1.0;
    }
}

void initializing_pressure(double P[N+2][N+2])
{
    for(int i = 0; i < N+2; i++)
    {
        for(int j = 0; j < N+2; j++)
        {
            P[i][j] = 0.0;
        }
    }
}

void intial_conditions_for_star_velocities(double u_star[N+2][N+2], double v_star[N+2][N+2])
{
    for(int i = 1; i <=N; i++)
    {
        for(int j = 1; j <=N; j++)
        {
            u_star[i][j] = 0.0;
            v_star[i][j] = 0.0;
        }
    }
    for(int i = 1; i<=N; i++)
    {
        // CORRECTED: Changed initial lid velocity for u_star from 2.0 to 1.0
        u_star[N+1][i] = 1.0; // top boundary ghost cells
        v_star[N+1][i] = 0.0;
        u_star[0][i] = 0.0; // Bottom boundary ghost cells
        v_star[0][i] = 0.0;
        u_star[i][0] = 0.0; // left boundary ghost cells
        v_star[i][0] = 0.0;
        u_star[i][N+1] = 0.0; // Right boundary ghost cells
        v_star[i][N+1] = 0.0;
    }
}

void initial_conditions_for_face_centre_velocities(double U[N][N+1], double V[N+1][N])
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N+1; j++)
        {
            U[i][j] = 0.0;
        }
    }
    for(int i = 0; i < N+1; i++)
    {
        for(int j = 0; j < N; j++)
        {
            V[i][j] = 0.0;
        }
    }
}

void compute_star_velocities_at_first_time_step(double Re, double dt, double dx, double dy, double aP, double aE, double aW, double aN, double aS, double u_star[N+2][N+2], double v_star[N+2][N+2], double u[N+2][N+2], double v[N+2][N+2], double U[N][N+1], double V[N+1][N])
{
    for(int i = 1; i <=N; i++)
    {
        for(int j = 1; j <=N; j++)
        {
            u_star[i][j] = pow(aP,-1)*(aE*u_star[i][j+1] + aW*u_star[i][j-1] + aN*u_star[i+1][j] + aS*u_star[i-1][j] + u[i][j]
                                        -(dt/dx)*(U[i-1][j]*0.5*(u[i][j+1] + u[i][j]) - U[i-1][j-1]*0.5*(u[i][j-1] + u[i][j]))
                                        -(dt/dy)*(V[i][j-1]*0.5*(u[i][j] + u[i+1][j]) - V[i-1][j-1]*0.5*(u[i][j] + u[i-1][j]))
                                        +(dt/(2*Re*pow(dx,2)))*(u[i][j+1] - 2*u[i][j] + u[i][j-1])
                                        +(dt/(2*Re*pow(dy,2)))*(u[i+1][j] - 2*u[i][j] + u[i-1][j]));

            v_star[i][j] = pow(aP,-1)*(aE*v_star[i][j+1] + aW*v_star[i][j-1] + aN*v_star[i+1][j] + aS*v_star[i-1][j] + v[i][j]
                                        -(dt/dx)*(U[i-1][j]*0.5*(v[i][j+1] + v[i][j]) - U[i-1][j-1]*0.5*(v[i][j-1] + v[i][j]))
                                        -(dt/dy)*(V[i][j-1]*0.5*(v[i+1][j] + v[i][j]) - V[i-1][j-1]*0.5*(v[i][j] + v[i-1][j]))
                                        +(dt/(2*Re*pow(dx,2)))*(v[i][j+1] - 2*v[i][j] + v[i][j-1])
                                        +(dt/(2*Re*pow(dy,2)))*(v[i+1][j] - 2*v[i][j] + v[i-1][j]));
        }
    }
}

void compute_star_velocities_after_first_time_step(double Re, double dt, double dx, double dy, double aP, double aE, double aW, double aN, double aS, double u_star[N+2][N+2], double v_star[N+2][N+2], double u[N+2][N+2], double u_prev[N+2][N+2], double v[N+2][N+2], double v_prev[N+2][N+2],double U[N][N+1], double U_prev[N][N+1], double V[N+1][N], double V_prev[N+1][N])
{
    for(int i = 1; i <=N; i++)
    {
        for(int j = 1; j <=N; j++)
        {
            u_star[i][j] = pow(aP,-1)*(aE*u_star[i][j+1] + aW*u_star[i][j-1] + aN*u_star[i+1][j] + aS*u_star[i-1][j] + u[i][j]
                                        -(dt/dx)*(U[i-1][j]*0.75*(u[i][j+1] + u[i][j]) - U[i-1][j-1]*0.75*(u[i][j-1] + u[i][j]))
                                        +(dt/dx)*(U_prev[i-1][j]*0.25*(u_prev[i][j+1] + u_prev[i][j]) - U_prev[i-1][j-1]*0.25*(u_prev[i][j-1] + u_prev[i][j]))
                                        -(dt/dy)*(V[i][j-1]*0.75*(u[i][j] + u[i+1][j]) - V[i-1][j-1]*0.75*(u[i][j] + u[i-1][j]))
                                        +(dt/dy)*(V_prev[i][j-1]*0.25*(u_prev[i][j] + u_prev[i+1][j]) - V_prev[i-1][j-1]*0.25*(u_prev[i][j] + u_prev[i-1][j]))
                                        +(dt/(2*Re*pow(dx,2)))*(u[i][j+1] - 2*u[i][j] + u[i][j-1])
                                        +(dt/(2*Re*pow(dy,2)))*(u[i+1][j] - 2*u[i][j] + u[i-1][j]));

            v_star[i][j] = pow(aP,-1)*(aE*v_star[i][j+1] + aW*v_star[i][j-1] + aN*v_star[i+1][j] + aS*v_star[i-1][j] + v[i][j]
                                        -(dt/dx)*(U[i-1][j]*0.75*(v[i][j+1] + v[i][j]) - U[i-1][j-1]*0.75*(v[i][j-1] + v[i][j]))
                                        +(dt/dx)*(U_prev[i-1][j]*0.25*(v_prev[i][j+1] + v_prev[i][j]) - U_prev[i-1][j-1]*0.25*(v_prev[i][j-1] + v_prev[i][j]))
                                        -(dt/dy)*(V[i][j-1]*0.75*(v[i+1][j] + v[i][j]) - V[i-1][j-1]*0.75*(v[i][j] + v[i-1][j]))
                                        +(dt/dy)*(V_prev[i][j-1]*0.25*(v_prev[i+1][j] + v_prev[i][j]) - V_prev[i-1][j-1]*0.25*(v_prev[i][j] + v_prev[i-1][j]))
                                        +(dt/(2*Re*pow(dx,2)))*(v[i][j+1] - 2*v[i][j] + v[i][j-1])
                                        +(dt/(2*Re*pow(dy,2)))*(v[i+1][j] - 2*v[i][j] + v[i-1][j]));
        }
    }
}


void update_star_velocities_of_ghost_cells(double u_star[N+2][N+2], double v_star[N+2][N+2])
{
    for(int i = 1; i <=N; i++)
    {
        u_star[i][0] = -u_star[i][1]; //Left boundary ghost cells
        u_star[i][N+1] = -u_star[i][N]; // Right biundary ghost cells
        v_star[i][0] = -v_star[i][1]; //Left boundary ghost cells
        v_star[i][N+1] = -v_star[i][N]; //Right biundary ghost cells
        u_star[0][i] = -u_star[1][i]; // Bottom
        u_star[N+1][i] = 2.0 - u_star[N][i]; // top (This is correct for lid velocity 1.0)
        v_star[0][i] = -v_star[1][i]; //bottom
        v_star[N+1][i] = -v_star[N][i]; //top
    }
}

void assign_starvelocities(double u_star[N+2][N+2], double v_star[N+2][N+2], double u_star_upd[N+2][N+2], double v_star_upd[N+2][N+2])
{
    for(int i = 1; i <=N; i++)
    {
        for(int j = 1; j <=N; j++)
        {
            u_star_upd[i][j] = u_star[i][j];
            v_star_upd[i][j] = v_star[i][j];
        }
    }
    for(int i = 1 ; i <= N; i++)
    {
        u_star_upd[i][0] = u_star[i][0]; // left boundary ghost cells
        u_star_upd[i][N+1] = u_star[i][N+1]; // right
        v_star_upd[i][0] = v_star[i][0]; // left
        v_star_upd[i][N+1] = v_star[i][N+1]; // right
        u_star_upd[0][i] = u_star[0][i]; // bottom
        u_star_upd[N+1][i] = u_star[N+1][i]; // top
        v_star_upd[0][i] = v_star[0][i]; // bottom
        v_star_upd[N+1][i] = v_star[N+1][i]; // top
    }
}

void update_previoustimestep_velocities(double u[N+2][N+2], double v[N+2][N+2], double u_prev[N+2][N+2], double v_prev[N+2][N+2])
{
    for(int i = 1; i <=N; i++)
    {
        for(int j = 1; j <=N; j++)
        {
            u_prev[i][j] = u[i][j];
            v_prev[i][j] = v[i][j];
        }
    }
    for(int i = 1 ; i <= N; i++)
    {
        u_prev[i][0] = u[i][0]; // left boundary ghost cells
        u_prev[i][N+1] = u[i][N+1]; // right
        v_prev[i][0] = v[i][0]; // left
        v_prev[i][N+1] = v[i][N+1]; // right
        u_prev[0][i] = u[0][i]; // bottom
        u_prev[N+1][i] = u[N+1][i]; // top
        v_prev[0][i] = v[0][i]; // bottom
        v_prev[N+1][i] = v[N+1][i]; // top
    }
}

void update_previoustimestep_face_centre_velocities(double U[N][N+1], double U_prev[N][N+1],double V_prev[N+1][N], double V[N+1][N])
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N+1; j++)
        {
            U_prev[i][j] = U[i][j];
        }
    }
    for(int i = 0; i < N+1; i++)
    {
        for(int j = 0; j < N; j++)
        {
            V_prev[i][j] = V[i][j];
        }
    }
}


void assign_Pressure_values(double P[N+2][N+2], double P_new[N+2][N+2])
{
    for(int i = 1; i <=N; i++)
    {
        for(int j = 1; j <=N; j++)
        {
            P_new[i][j] = P[i][j];
            P_new[i][j] = P[i][j]; // This line is redundant but kept as per original
        }
    }
    for(int i = 1 ; i <= N; i++)
    {
       P_new[i][0] = P[i][0]; // left boundary ghost cells
       P_new[i][N+1] = P[i][N+1]; // right
       P_new[i][0] = P[i][0]; // left (This line is redundant but kept as per original)
       P_new[i][N+1] = P[i][N+1]; // right (This line is redundant but kept as per original)
       P_new[0][i] = P[0][i]; // bottom
       P_new[N+1][i] = P[N+1][i]; // top
       P_new[0][i] = P[0][i]; // bottom (This line is redundant but kept as per original)
       P_new[N+1][i] = P[N+1][i]; // top (This line is redundant but kept as per original)
    }
}

void solve_pressure_poisson_equation(double dx, double dt, double P[N+2][N+2], double u_star[N+2][N+2], double v_star[N+2][N+2])
{
    for(int i = 1; i <=N; i++)
    {
        for(int j = 1; j <=N; j++)
        {
            P[i][j] = 0.25*(P[i][j+1] + P[i][j-1] + P[i+1][j] + P[i-1][j] - (dx/dt)*0.5*(u_star[i][j+1] - u_star[i][j-1]) - (dx/dt)*0.5*(v_star[i+1][j] - v_star[i-1][j]));

        }
    }
}

void update_pressure_for_ghost_cells(double P[N+2][N+2])
{
    for(int i = 1; i <=N; i++)
    {
        P[0][i] = P[1][i];
        P[N+1][i] = P[N][i];
        P[i][0] = P[i][1];
        P[i][N+1] = P[i][N];
    }
}

void compute_velocities_at_cell_centres(double dt, double dx, double dy, double u[N+2][N+2], double v[N+2][N+2], double u_star[N+2][N+2], double v_star[N+2][N+2], double P[N+2][N+2])
{
    for(int i = 1; i <=N; i++)
    {
        for(int j = 1; j <=N; j++)
        {
            u[i][j] = u_star[i][j] - (dt/dx)*0.5*(P[i][j+1] - P[i][j-1]);
            v[i][j] = v_star[i][j] - (dt/dx)*0.5*(P[i+1][j] - P[i-1][j]);
        }
    }
}

void update_velocities_of_ghostcells(double u[N+2][N+2], double v[N+2][N+2])
{
    for(int i = 1; i <=N; i++)
    {
        u[i][0] = -u[i][1];
        u[i][N+1] = -u[i][N];
        v[i][0] = -v[i][1];
        v[i][N+1] = -v[i][N];
        u[0][i] = -u[1][i];
        u[N+1][i] = 2.0 - u[N][i]; // top (This is correct for lid velocity 1.0)
        v[0][i] = -v[1][i];
        v[N+1][i] = -v[N][i];
    }
}

// void compute_velocities_at_face_centres(double dt, double dx, double dy, double U[N][N+1], double V[N+1][N], double u_star[N+2][N+2], double v_star[N+2][N+2], double P[N+2][N+2])
// {
//     for(int i = 0; i < N; i++)
//     {
//         for(int j = 0; j <=N; j++)
//         {
//             U[i][j] = 0.5*(u_star[i+1][j+1] + u_star[i][j+1]) - (dt/dx)*(P[i+1][j+1] - P[i][j+1]);
//         }
//     }
//     for(int i = 0; i<=N; i++)
//     {
//         for(int j = 0; j < N; j++)
//         {
//             V[i][j] = 0.5*(v_star[i+1][j+1] + v_star[i+1][j]) - (dt/dy)*(P[i+1][j+1] - P[i+1][j]);
//         }
//     }
// }

void compute_velocities_at_face_centres(double dt, double dx, double dy, double U[N][N+1], double V[N+1][N], double u[N+2][N+2], double v[N+2][N+2])
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j <=N; j++)
        {
            U[i][j] = (u[i][j+1] + u[i+1][j+1])/2.0; //0.5*(u_star[i+1][j+1] + u_star[i+1][j]) - (dt/dx)*(P[i+1][j+1] - P[i+1][j]);
        }
    }
    for(int i = 0; i<=N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            V[i][j] = (v[i+1][j] + v[i+1][j+1])/2.0; //0.5*(v_star[i+1][j+1] + v_star[i][j+1]) - (dt/dy)*(P[i+1][j+1] - P[i][j+1]);
        }
    }
}

double infinity_norm(double u[N+2][N+2], double v[N+2][N+2], double u_upd[N+2][N+2], double v_upd[N+2][N+2])
{
    double max_value = 0.0;
    for(int i = 1; i<=N; i++)
    {
        for(int j = 1; j<=N; j++)
        {
            if(fabs(u[i][j] - u_upd[i][j]) > max_value)
            {
                max_value = fabs(u[i][j] - u_upd[i][j]);
            }
        }
    }
    for(int i = 1; i<=N; i++)
    {
        for(int j = 1; j<=N; j++)
        {
            if(fabs(v[i][j] - v_upd[i][j]) > max_value)
            {
                max_value = fabs(v[i][j] - v_upd[i][j]);
            }
        }
    }
    return max_value;
}

double infinity_norm_for_pressure(double u[N+2][N+2], double u_upd[N+2][N+2])
{
    double max_value = 0.0;
    for(int i = 1; i<=N; i++)
    {
        for(int j = 1; j<=N; j++)
        {
            if(fabs(u[i][j] - u_upd[i][j]) > max_value)
            {
                max_value = fabs(u[i][j] - u_upd[i][j]);
            }
        }
    }
    return max_value;
}


int main()
{
    //int N = 128;
    double dx = 1/(double)N;
    double dy = 1/(double)N;
    int Re = 100;
    double CFL = 0.3;
    double dt = CFL*dx;
    double sim_time = 0.0;
    double u[N+2][N+2], v[N+2][N+2], u_star[N+2][N+2], v_star[N+2][N+2], u_star_upd[N+2][N+2], v_star_upd[N+2][N+2], u_prev[N+2][N+2], v_prev[N+2][N+2];
    double U[N][N+1], U_prev[N][N+1], U_star[N][N+1], V[N+1][N], V_prev[N+1][N], V_star[N+1][N];
    double P[N+2][N+2], P_new[N+2][N+2];
    double aP, aE, aW, aN, aS;
    double norm = 1;
    double norm_P = 1;
    int iter = 0;
    aP = 1 + (2*dt)/(Re*pow(dx,2));
    aE = dt/(2*Re*pow(dx,2));
    aW = aE;
    aN = aE;
    aS = aE;
    intial_conditions_for_cell_center_velocities_and_Pressure(u, v, P);
    initial_conditions_for_face_centre_velocities(U, V);
    intial_conditions_for_star_velocities(u_star, v_star);
    do
    {
        assign_starvelocities(u_star, v_star, u_star_upd, v_star_upd);
        compute_star_velocities_at_first_time_step(Re,dt,dx,dy,aP,aE,aW,aN,aS,u_star,v_star,u,v,U,V);
        //update_star_velocities_of_ghost_cells(u_star, v_star); // This line was commented in your original, keeping it commented.
        norm = infinity_norm(u_star, v_star, u_star_upd, v_star_upd);
    } while (norm > epsilon);
    update_star_velocities_of_ghost_cells(u_star, v_star);
    do
    {
        assign_Pressure_values(P, P_new);
        solve_pressure_poisson_equation(dx, dt, P, u_star, v_star);
        //update_pressure_for_ghost_cells(P); // This line was commented in your original, keeping it commented.
        norm_P = infinity_norm_for_pressure(P, P_new);
        iter++;
       //cout << "iter = " << iter << "   norm = " << norm_P << endl;
    } while (norm_P > epsilon);
    update_pressure_for_ghost_cells(P);
    compute_velocities_at_cell_centres(dt, dx, dy, u, v, u_star, v_star, P);
    update_velocities_of_ghostcells(u,v);
    //compute_velocities_at_face_centres(dt, dx, dy, U, V, u_star, v_star, P);
    compute_velocities_at_face_centres(dt, dx, dy, U, V, u, v);

    intial_conditions_for_cell_center_velocities_and_Pressure(u_prev, v_prev, P);
    initial_conditions_for_face_centre_velocities(U_prev, V_prev);
    //intial_conditions_for_star_velocities(u_star, v_star); // This line was commented in your original, keeping it commented.
    do
    {
        iter = 0;
        intial_conditions_for_star_velocities(u_star, v_star);
        //initializing_pressure(P); // This line was commented in your original, keeping it commented.
        assign_Pressure_values(P_new, P); // Use previous pressure as initial guess
        do
        {
            assign_starvelocities(u_star, v_star, u_star_upd, v_star_upd);
            compute_star_velocities_after_first_time_step(Re,dt,dx,dy,aP,aE,aW,aN,aS,u_star,v_star,u,u_prev,v,v_prev,U,U_prev,V,V_prev);
            norm = infinity_norm(u_star, v_star, u_star_upd, v_star_upd);

        } while (norm > epsilon);
        update_star_velocities_of_ghost_cells(u_star, v_star);
        update_previoustimestep_velocities(u,v,u_prev,v_prev);
        update_previoustimestep_face_centre_velocities(U, U_prev, V_prev, V);

        //update_pressure_for_ghost_cells(P); // This line was commented in your original, keeping it commented.
        do
        {
            assign_Pressure_values(P, P_new);
            solve_pressure_poisson_equation(dx, dt, P, u_star, v_star);
            norm_P = infinity_norm_for_pressure(P, P_new);
            iter++;

        } while (norm_P > epsilon);

        update_pressure_for_ghost_cells(P);
        compute_velocities_at_cell_centres(dt, dx, dy, u, v, u_star, v_star, P);
        update_velocities_of_ghostcells(u,v);
        //compute_velocities_at_face_centres(dt, dx, dy, U, V, u_star, v_star, P);
        compute_velocities_at_face_centres(dt, dx, dy, U, V, u, v);

        sim_time = sim_time + dt;
        cout << "\n simulation time = " << sim_time << endl;
        cout << "No.of iterations taken for pressure poisson equation to converge for each time step = " << iter << endl;
        cout << "Norm for Pressure = " << norm_P << endl;


    } while (sim_time < 7.0);


    for(int i = N; i>=1; i--)
    {
        cout << u[i][64] << endl;
    }

    cout << "\n simulation time = " << sim_time << endl;
    return 0;
}