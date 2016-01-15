#ifndef ENVELOPE_H
#define ENVELOPE_H
class Envelope
{
public:
    double xmin;
    double ymin;
    double xmax;
    double ymax;
    Envelope(double xmin,double ymin,double xmax,double ymax):xmin(xmin),ymin(ymin),xmax(xmax),ymax(ymax){}
};

#endif // ENVELOPE_H
