#include <stdio.h>
#include <math.h>
#include <locale.h>
#include <omp.h>
#include <stdlib.h>
typedef struct {
    double mass;
    double radius;
    double position[3];
    double velocity[3];
} Body;

double distance(Body body1, Body body2) {
    double dx = body1.position[0] - body2.position[0];
    double dy = body1.position[1] - body2.position[1];
    double dz = body1.position[2] - body2.position[2];
    return sqrt(dx * dx + dy * dy + dz * dz);
}

void calculateGravitationalForce(Body* body1, Body* body2, double G) {
    double r = distance(*body1, *body2);
    double F = G * body1->mass * body2->mass / (r * r);

    for (int i = 0; i < 3; i++) {
        double dx = body2->position[i] - body1->position[i];
        body1->velocity[i] += F * dx / (r * body1->mass);
        //body2->velocity[i] -= F * dx / (r * body2->mass);
    }
}

void updateBody(Body* body, double dt) {
    for (int i = 0; i < 3; i++) {
        body->position[i] += body->velocity[i] * dt;
    }
}

void handleElasticCollision(Body* body1, Body* body2) {
    double r = distance(*body1, *body2);
    if (r < body1->radius + body2->radius) {
        // Вычисляем нормализованный вектор разницы позиций
        double n[3];
        for (int i = 0; i < 3; i++) {
            n[i] = (body2->position[i] - body1->position[i]) / r;
        }
            
        // Вычисляем относительную скорость
        double v_rel[3];
        for (int i = 0; i < 3; i++) {
            v_rel[i] = body2->velocity[i] - body1->velocity[i];
        }

        // Вычисляем скорости после соударения (упругое соударение)
        double impulse = 2.0 * body1->mass * body2->mass / (body1->mass + body2->mass) * (v_rel[0] * n[0] + v_rel[1] * n[1] + v_rel[2] * n[2]);
#pragma omp critical
        {
        for (int i = 0; i < 3; i++) {
            body1->velocity[i] += impulse / body1->mass * n[i];
            body2->velocity[i] -= impulse / body2->mass * n[i];
        }

        // Избегаем "слипания" тел
        double overlap = (body1->radius + body2->radius - r) / 2.0;
        for (int i = 0; i < 3; i++) {
            body1->position[i] -= overlap * n[i];
            body2->position[i] += overlap * n[i];
        }
        }
    }
}

void movementBodies(Body* bodies,double G, double dt, double totalTime, int numBodies){
for (double t = 0; t < totalTime; t += dt) {
        // Вычисление гравитационных сил и обновление положения
#pragma omp parallel for
            for (int i = 0; i < numBodies; i++) {
                for (int j = 0; j < numBodies; j++) {
                    if (i != j) {
                        calculateGravitationalForce(&bodies[i], &bodies[j], G);
                    }
                }
            }
        // Обновление значений массива тел
//#pragma omp parallel for
        for (int i = 0; i < numBodies; i++) {
            updateBody(&bodies[i], dt);
        }

        // Обработка упругих соударений
#pragma omp parallel for
        for (int i = 0; i < numBodies; i++) {
            for (int j = i + 1; j < numBodies; j++) {
                handleElasticCollision(&bodies[i], &bodies[j]);
            }
        }
    }
}
Body initValues(double mass, double radius, double x, double y, double z, double vx, double vy, double vz){
    Body body;
    body.mass = mass;
    body.radius  = radius;
	body.position[0] = x;
	body.position[1] = y;
	body.position[2] = z;
	body.velocity[0] = vx;
	body.velocity[1] = vy;
	body.velocity[2] = vz;
    return body;
}
void printBody(Body body){
    printf("Масса, радиус, положение и скорость тела :(%.1lf), (%.2lf) (%lf, %lf, %lf), (%lf, %lf, %lf)\n", body.mass, body.radius,
            body.position[0], body.position[1], body.position[2],
            body.velocity[0], body.velocity[1], body.velocity[2]);
}
int main() {
    //printf("%d\n",omp_get_max_threads());
    //omp_set_num_threads(omp_get_max_threads());
    setlocale(LC_ALL, "Rus");
    int numBodies = 2000;  // Количество тел
    Body* bodies;
    double G = 6.67e-11; // гравитационная постоянная
    double totalTime = 3.5;
    double dt = 0.01; // шаг времени
    int printResult = 0;
	while(1){
	
	int tmp;
	printf("G=%g\tdt=%f\tN_body=%d\t printResult=%d\n",G,dt,numBodies,printResult);
	printf("you need change setting? (0 - no; 1 - yes; 10 -exit; 11 - testing mode)\n");
	scanf("%d",&tmp);
	if(tmp == 1) {
	printf("Enter G, dt, numBodies, printResult(1/0) (Example: 10 0.01 1000 0)\n");
	scanf("%lf %lf %d %d",&G,&dt,&numBodies,&printResult);
	}
	if(tmp == 10) break;
    printf("Enter totalTime (Example: 1.5)\n");
    scanf("%lf", &totalTime);

    // Инициализация параметров тел
    if(tmp == 11){
    numBodies = 3;
    //G = 6.67e-11;
    printResult = 1;
    dt = 0.01;
    bodies = (Body*)malloc(sizeof(Body)*numBodies);
    bodies[0] = initValues(10.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    bodies[1] = initValues(5.0, 0.5, 5.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    bodies[2] = initValues(500.0, 0.5, -5.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    else{
        bodies = (Body*)malloc(sizeof(Body)*numBodies);
        for (int i = 0; i < numBodies; i++) {
            bodies[i] = initValues(10.0 + rand()%100, (double)(1 + rand() % 4), 
            (double)(rand()%100), (double)(rand() % 100), (double)i * 10, 
            10-(double)(rand() % 20), 10-(double)(rand() % 20), 10-(double)(rand() % 20));
        }
    }
    for(int i = 0; i < numBodies && printResult; i++){
        if(!i) printf("Начальные значения:\n\n");
        printf("%d. ",i);
        printBody(bodies[i]);
    }
    double start = omp_get_wtime();

    movementBodies(bodies, G, dt, totalTime, numBodies);

    double end = omp_get_wtime();
    for (int i = 0; i < numBodies && printResult; i++) {
        if(!i) printf("\nКонечные значения:\n\n");
        printf("%d. ",i);
        printBody(bodies[i]);
        //printf("Положение и скорость тела %d: (%f, %f, %f), (%f, %f, %f)\n", i + 1,
          //  bodies[i].position[0], bodies[i].position[1], bodies[i].position[2],
           // bodies[i].velocity[0], bodies[i].velocity[1], bodies[i].velocity[2]);
    }
    printf("\nTIME: %f (sec)\n\n\n", end-start);
    free(bodies);
    }
    return 0;
}

