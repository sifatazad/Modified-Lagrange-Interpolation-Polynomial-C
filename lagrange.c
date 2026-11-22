#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define Abs(a) ((a)>=0?(a):-(a))

struct data {
    double x;
    double y;
};

static int compare(const void *a, const void *b) {
    const struct data *p= (const struct data*) a;
    const struct data *q= (const struct data*) b;

    if (p->x > q->x)
        return 1;
    if (p->x < q->x)
        return -1;
    return 0;
}

static int total_point(FILE *fp) {
    int n=0;
    char line[500];
    while (fgets(line,500,fp)!=NULL) {
        n++;
    }
    rewind(fp);
    return n;
}

static double lagrange(FILE *fp, const double x_r) {
    const int n= total_point(fp);

    if (n<2) {
        printf("Error! At least two points are required.\n");
        return NAN;
    }

    double *x= malloc(sizeof(double)*n);
    if (!x) {
        return NAN;
    }

    double *y= malloc(sizeof(double)*n);
    if (!y) {
        free(x);
        return NAN;
    }

    for (int i=0; i<n; i++) {
        fscanf(fp,"%lf\t%lf",&x[i],&y[i]);
    }
    rewind(fp);

    struct data *point= malloc(sizeof(struct data)*n);
    if (!point) {
        free(x);
        free(y);
        return NAN;
    }

    for (int i=0; i<n; i++) {
        (point+i)->x= x[i];
        (point+i)->y= y[i];
    }

    qsort(point,n,sizeof(struct data),compare);

    for (int i=0; i<n; i++) {
        x[i]=(point+i)->x;
        y[i]=(point+i)->y;
    }

    if (x_r<x[0] || x_r>x[n-1]) {
        printf("Error! Extrapolation not possible!\n");
        free(x);
        free(y);
        free(point);
        return NAN;
    }

    double s=0;
    for (int i=0; i<n; i++) {
        double mul=1;
        for (int j=0; j<n; j++) {
            if (i==j) {
                continue;
            }
            mul*=(x_r-x[j])/(x[i]-x[j]);
        }
        mul*=y[i];
        s+=mul;
    }

    free(x);
    free(y);
    free(point);
    return s;
}

static double lagrange_without_file(double *x, double *y, const double x_r, const int n) {
    if (n<2) {
        printf("Error! At least two points are required.\n");
        return NAN;
    }

    struct data *point= malloc(sizeof(struct data)*n);
    if (!point) {
        return NAN;
    }

    for (int i=0; i<n; i++) {
        (point+i)->x= x[i];
        (point+i)->y= y[i];
    }

    qsort(point,n,sizeof(struct data),compare);

    for (int i=0; i<n; i++) {
        x[i]=(point+i)->x;
        y[i]=(point+i)->y;
    }

    if (x_r<x[0] || x_r>x[n-1]) {
        printf("Error! Extrapolation not possible!\n");
        free(point);
        return NAN;
    }

    double s=0;
    for (int i=0; i<n; i++) {
        double mul=1;
        for (int j=0; j<n; j++) {
            if (i==j) {
                continue;
            }
            mul*=(x_r-x[j])/(x[i]-x[j]);
        }
        mul*=y[i];
        s+=mul;
    }

    free(point);
    return s;
}

static double m_lagrange(FILE *fp, const double x_r) {
    const int n= total_point(fp);

    if (n<2) {
        printf("Error! At least two points are required.\n");
        return NAN;
    }

    if (n>=6) {
        double *x= malloc(sizeof(double)*n);
        if (!x) {
            return NAN;
        }

        double *y= malloc(sizeof(double)*n);
        if (!y) {
            free(x);
            return NAN;
        }

        for (int i=0; i<n; i++) {
            fscanf(fp,"%lf\t%lf",&x[i],&y[i]);
        }
        rewind(fp);

        struct data *point= malloc(sizeof(struct data)*n);
        if (!point) {
            free(x);
            free(y);
            return NAN;
        }

        for (int i=0; i<n; i++) {
            (point+i)->x= x[i];
            (point+i)->y= y[i];
        }

        qsort(point,n,sizeof(struct data),compare);

        for (int i=0; i<n; i++) {
            x[i]=(point+i)->x;
            y[i]=(point+i)->y;
        }

        if (x_r<x[0] || x_r>x[n-1]) {
            printf("Error! Extrapolation not possible!\n");
            free(x);
            free(y);
            free(point);
            return NAN;
        }

        int interval= 0;
        for (int i=0; i<n; i++) {
            if (x[i]>x_r) {
                interval=i-1;
                break;
            }
            else if (x_r==x[i]) {
                double val= y[i];
                free(x);
                free(y);
                free(point);
                return val;
            }
        }

        const int left= interval;
        const int right= n-left-2;

        double *new_x= malloc(sizeof(double)*6);
        if (!new_x) {
            free(x);
            free(y);
            free(point);
            return NAN;
        }

        double *new_y= malloc(sizeof(double)*6);
        if (!new_y) {
            free(x);
            free(y);
            free(new_x);
            free(point);
            return NAN;
        }

        if (left>=2 && right>=2) {
            for (int i=interval-2; i<=interval+3; i++) {
                new_x[i-(interval-2)]= x[i];
                new_y[i-(interval-2)]= y[i];
            }
        }
        else if (left<2) {
            for (int i=0; i<6; i++) {
                new_x[i]= x[i];
                new_y[i]= y[i];
            }
        }
        else if (right<2) {
            for (int i=n-6; i<n; i++) {
                new_x[i-n+6]= x[i];
                new_y[i-n+6]= y[i];
            }
        }

        double val= lagrange_without_file(new_x, new_y, x_r, 6);

        free(x);
        free(y);
        free(new_x);
        free(new_y);
        free(point);
        return val;
    }
    else {
        const double val= lagrange(fp,x_r);
        return val;
    }
}

static double m_lagrange_without_file(double *x, double *y, const double x_r, const int n) {
    if (n<2) {
        printf("Error! At least two points are required.\n");
        return NAN;
    }

    if (n>=6) {
        struct data *point= malloc(sizeof(struct data)*n);
        if (!point) {
            return NAN;
        }

        for (int i=0; i<n; i++) {
            (point+i)->x= x[i];
            (point+i)->y= y[i];
        }

        qsort(point,n,sizeof(struct data),compare);

        for (int i=0; i<n; i++) {
            x[i]=(point+i)->x;
            y[i]=(point+i)->y;
        }

        if (x_r<x[0] || x_r>x[n-1]) {
            printf("Error! Extrapolation not possible!\n");
            free(point);
            return NAN;
        }

        int interval= 0;
        for (int i=0; i<n; i++) {
            if (x[i]>x_r) {
                interval=i-1;
                break;
            }
            else if (x_r==x[i]) {
                double val= y[i];
                free(point);
                return val;
            }
        }

        const int left= interval;
        const int right= n-left-2;

        double *new_x= malloc(sizeof(double)*6);
        if (!new_x) {
            free(point);
            return NAN;
        }

        double *new_y= malloc(sizeof(double)*6);
        if (!new_y) {
            free(new_x);
            free(point);
            return NAN;
        }

        if (left>=2 && right>=2) {
            for (int i=interval-2; i<=interval+3; i++) {
                new_x[i-(interval-2)]= x[i];
                new_y[i-(interval-2)]= y[i];
            }
        }
        else if (left<2) {
            for (int i=0; i<6; i++) {
                new_x[i]= x[i];
                new_y[i]= y[i];
            }
        }
        else if (right<2) {
            for (int i=n-6; i<n; i++) {
                new_x[i-n+6]= x[i];
                new_y[i-n+6]= y[i];
            }
        }

        double val= lagrange_without_file(new_x, new_y, x_r, 6);

        free(new_x);
        free(new_y);
        free(point);
        return val;
    }
    else {
        const double val= lagrange_without_file(x,y,x_r,n);
        return val;
    }
}

int main() {
    FILE *fp= fopen("C:/Users/Rifat/Downloads/50_data_points.txt","r");
    if (!fp) {
        return 1;
    }

    int n=total_point(fp);
    printf("%.6lf\n",m_lagrange(fp,9.10));

    fclose(fp);
}