#include <mystdlib.h>

#include <myadt.hpp>
#include <meshing.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>

#include "fieldlines.hpp"

namespace netgen
{
  // Along a field line the next point is almost always inside the element we
  // just left, or one next to it. Test that element first: a hit costs one
  // PointContainedIn3DElement, a miss falls back to the search tree (which
  // also allocates). Mesh::SetPointSearchStartElement looks like it does this
  // already, but ps_startelement is never read by the search.
  inline int GetVolElementHinted(const Mesh& mesh, const Point<3>& p,
                                 double* lami, int hint);

  inline int GetVolElement(const Mesh& mesh, const Point<3>& p,
                           double* lami)
  {
    if(mesh.GetDimension() == 3)
      {
        auto ei = mesh.GetElementOfPoint(p, lami, true);
        if(!ei.IsValid())
          return -1;
        return ei.Nr0();
      }
    else
      {
        auto ei = mesh.GetSurfaceElementOfPoint(p, lami, true);
        if(!ei.IsValid())
          return -1;
        return ei.Nr0();
      }
  }

  inline int GetVolElementHinted(const Mesh& mesh, const Point<3>& p,
                                 double* lami, int hint)
  {
    if(hint >= 0 && mesh.GetDimension() == 3)
      {
        ElementIndex ei = ElementIndex::FromNr0(hint);
        if(mesh.PointContainedIn3DElement(p, lami, ei))
          return hint;
      }
    return GetVolElement(mesh, p, lami);
  }

  RKStepper :: ~RKStepper() 
  {
    delete a;
  }
    
  RKStepper :: RKStepper(int type) : a(NULL), tolerance(1e100)
  {
    notrestarted = 0;

    if (type == 0) // explicit Euler
      {
        c.SetSize(1); c[0] = 0;
        b.SetSize(1); b[0] = 1;
        steps = order = 1;
      }
    else if (type == 1) // Euler-Cauchy
      {
        c.SetSize(2); c[0] = 0; c[1] = 0.5;
        b.SetSize(2); b[0] = 0; b[1] = 1;
        Array<int> size(2);
        size[0] = 0; size[1] = 1;
        a = new TABLE<double>(size);
        a->Set(2,1,0.5);  // Set, Get: 1-based!
        steps = order = 2;
      }
    else if (type == 2) // Simpson
      {
        c.SetSize(3); c[0] = 0; c[1] = 1; c[2] = 0.5;
        b.SetSize(3); b[0] = b[1] = 1./6.; b[2] = 2./3.;
        Array<int> size(3);
        size[0] = 0; size[1] = 1; size[2] = 2;
        a = new TABLE<double>(size);
        a->Set(2,1,1);
        a->Set(3,1,0.25); a->Set(3,2,0.25); 
        steps = order = 3;
      }
    else if (type == 3) // classical Runge-Kutta
      {
        c.SetSize(4); c[0] = 0; c[1] = c[2] = 0.5; c[3] = 1;
        b.SetSize(4); b[0] = b[3] = 1./6.; b[1] = b[2] = 1./3.;
        Array<int> size(4);
        size[0] = 0; size[1] = 1; size[2] = 2; size[3] = 3;
        a = new TABLE<double>(size);
        a->Set(2,1,0.5);
        a->Set(3,1,0); a->Set(3,2,0.5); 
        a->Set(4,1,0); a->Set(4,2,0); a->Set(4,3,1); 
        steps = order = 4;
      }
    
    K.SetSize(steps);
  }

  void RKStepper :: StartNextValCalc(const Point<3> & astartval, const double astartt, const double ah, const bool aadaptive)
  {
    //cout << "Starting RK-Step with h=" << ah << endl;

    stepcount = 0;
    h = ah;
    startt = astartt;
    startval = astartval;
    adaptive = aadaptive;
    adrun = 0;
  }

  bool RKStepper :: GetNextData(Point<3> & val, double & t, double & ah)
  {
    bool finished = false;
    
    if(stepcount <= steps && stepcount>0)
      {
        t = startt + c[stepcount-1]*h;
        val = startval;
        for(int i=0; i<stepcount-1; i++)
          val += h * a->Get(stepcount,i+1) * K[i];
      }
    
    
    if(stepcount == steps)
      {
        val = startval;
        for(int i=0; i<steps; i++)
          val += h * b[i] * K[i];
        
        if(adaptive)
          {
            if(adrun == 0)
              {
                stepcount = 0;
                h *= 0.5;
                adrun = 1;
                valh = val;
              }
            else if (adrun == 1)
              {
                stepcount = 0;
                startval_bak = startval;
                startval = val;
                startt_bak = startt;
                startt += h;//0.5*h;
                adrun = 2;
              }
            else if (adrun == 2)
              {
                Point<3> valh2 = val;
                val = valh2 + 1./(pow(2.,order)-1.) * (valh2 - valh);
                auto errvec = val - valh;
                
                double err = errvec.Length();
                
                double fac = 0.7 * pow(tolerance/err,1./(order+1.));
                if(fac > 1.3) fac = 1.3;
                
                if(fac < 1 || notrestarted >= 2)
                  ah = 2.*h * fac;
                
                if(err < tolerance) 
                  {
                    finished = true;
                    notrestarted++;
                    //(*testout) << "finished RK-Step, new h=" << ah << " tolerance " << tolerance << " err " << err << endl;
                  }
                else
                  {
                    //ah *= 0.9;
                    notrestarted = 0;
                    //(*testout) << "restarting h " << 2.*h << " ah " << ah << " tolerance " << tolerance << " err " << err << endl;
                    StartNextValCalc(startval_bak,startt_bak, ah, adaptive);
                  }
              }
          }
        else 
          {
            t = startt + h;
            finished = true;
          }
        
      }
    
    if(stepcount == 0)
      {
        t = startt + c[stepcount]*h;
        val = startval;
        for(int i=0; i<stepcount; i++)
          val += h * a->Get(stepcount,i) * K[i];
      }
    
    return finished;
  }


  bool RKStepper :: FeedNextF(const Vec<3> & f)
  {
    K[stepcount] = f;
    stepcount++;
    return true;
  }
  


  void FieldLineCalc :: GenerateFieldLines(Array<Point<3>> & potential_startpoints, const int numlines)
  {
    pstart.SetSize0();
    pend.SetSize0();
    values.SetSize0();

    const int nstart = potential_startpoints.Size();
    double crit = 1.0;

    if(randomized)
      {
        double sum = 0;
        double lami[3];
        Vec<3> v;
        int hint = -1;

        for(int i=0; i<nstart; i++)
          {
            int elnr = GetVolElementHinted(mesh, potential_startpoints[i], lami, hint);
            if (elnr == -1)
              continue;
            hint = elnr;

            func(elnr, lami, v);
            sum += v.Length();
          }

        crit = sum/double(numlines);
      }

    // Draw every threshold up front, in the original order, so the result
    // does not depend on how the work is distributed over threads.
    Array<double> crit_values(nstart);
    for(int i=0; i<nstart; i++)
      crit_values[i] = randomized ? (double(rand())/RAND_MAX)*crit : critical_value;

    // Lines are independent, so trace them in parallel. Work is done in
    // chunks and merged in index order: the output is identical to the
    // sequential version, and we stop as soon as numlines lines are usable
    // instead of tracing every candidate.
    const int nthreads = max(1, ngcore::TaskManager::GetNumThreads());
    const int chunk = max(1, 4*nthreads);

    Array<Array<Point<3>>> chunk_pstart(chunk), chunk_pend(chunk);
    Array<Array<double>> chunk_values(chunk);

    int calculated = 0;

    for(int base = 0; base < nstart && calculated < numlines; base += chunk)
      {
        const int n = min2(chunk, nstart-base);

        ngcore::ParallelFor(n, [&](size_t k)
          {
            chunk_pstart[k].SetSize0();
            chunk_pend[k].SetSize0();
            chunk_values[k].SetSize0();

            RKStepper stepper(rk_type);
            stepper.SetTolerance(tolerance);

            Array<Point<3>> line_points;
            Array<double> line_values;
            Array<bool> drawelems;
            Array<int> dirstart;

            Calc(potential_startpoints[base+k], line_points, line_values,
                 drawelems, dirstart, stepper, crit_values[base+k]);

            for(int j=1; j<dirstart.Size(); j++)
              for(int l=dirstart[j-1]; l<dirstart[j]-1; l++)
                {
                  if(!drawelems[l] || !drawelems[l+1]) continue;

                  chunk_pstart[k].Append(line_points[l]);
                  chunk_pend[k].Append(line_points[l+1]);
                  chunk_values[k].Append( 0.5*(line_values[l]+line_values[l+1]) );
                }
          });

        for(int k = 0; k < n && calculated < numlines; k++)
          {
            if(chunk_pstart[k].Size() == 0) continue;

            for(auto l : Range(chunk_pstart[k]))
              {
                pstart.Append(chunk_pstart[k][l]);
                pend.Append(chunk_pend[k][l]);
                values.Append(chunk_values[k][l]);
              }
            calculated++;
          }
      }
  }



  FieldLineCalc :: FieldLineCalc(const Mesh & amesh, const VectorFunction & afunc,
                                 const double rel_length, const int amaxpoints, 
                                 const double rel_thickness, const double rel_tolerance, const int rk_type, const int adirection) :
    mesh(amesh), func(afunc), stepper(rk_type)
  {
    mesh.GetBox (pmin, pmax);
    rad = 0.5 * Dist (pmin, pmax);
    

    maxlength = (rel_length > 0) ? rel_length : 0.5;
    maxlength *= 2.*rad;

    thickness = (rel_thickness > 0) ? rel_thickness : 0.0015;
    thickness *= 2.*rad;
    
    double auxtolerance = (rel_tolerance > 0) ? rel_tolerance : 1.5e-3;
    auxtolerance *= 2.*rad;

    stepper.SetTolerance(auxtolerance);
    this->rk_type = rk_type;
    this->tolerance = auxtolerance;

    direction = adirection;
    
    
    maxpoints = amaxpoints;

    if(direction == 0)
      {
        maxlength *= 0.5;
        maxpoints /= 2;
      }
    

    critical_value = -1;

    randomized = false;
    
  }
  

  FieldLineCalc :: ~FieldLineCalc() {;}

  
  void FieldLineCalc :: Calc(const Point<3> & startpoint, Array<Point<3>> & points, Array<double> & vals, Array<bool> & drawelems, Array<int> & dirstart)
  {
    Calc(startpoint, points, vals, drawelems, dirstart, stepper, critical_value);
  }

  // Re-entrant: the stepper and the critical value are arguments, and no
  // mesh state is written, so several lines can be traced at once.
  void FieldLineCalc :: Calc(const Point<3> & startpoint, Array<Point<3>> & points, Array<double> & vals,
                             Array<bool> & drawelems, Array<int> & dirstart,
                             RKStepper & stepper, double crit_value) const
  {
    Vec<3> v = 0.0;
    double startlami[3] = {0.0, 0.0, 0.0};
    
    points.SetSize(0);
    vals.SetSize(0);
    drawelems.SetSize(0);

    dirstart.SetSize(0);
    dirstart.Append(0);

    int startelnr = GetVolElement(mesh, startpoint,startlami);
    if (startelnr == -1)
      return;
      

    Vec<3> startv;
    bool startdraw = func(startelnr, startlami, startv);

    double startval = startv.Length();

    if(crit_value > 0 && fabs(startval) < crit_value)
      return;

    //cout << "p = " << startpoint << "; elnr = " << startelnr << endl;


      
    for(int dir = 1; dir >= -1; dir -= 2)
      {
        if(dir*direction < 0) continue;
          
        points.Append(startpoint);
        vals.Append(startval);
        drawelems.Append(startdraw);
          
        double h = 0.001*rad/startval; // otherwise no nice lines; should be made accessible from outside
        
        v = startv;
        if(dir == -1) v *= -1.;

        int elnr = startelnr;
        double lami[3] = { startlami[0], startlami[1], startlami[2]}; 
          

        for(double length = 0; length < maxlength; length += h*vals.Last())
          {
            if(v.Length() < 1e-12*rad)
              {
                break;
              }

            double dummyt{0};
            stepper.StartNextValCalc(points.Last(),dummyt,h,true);
            stepper.FeedNextF(v);
            bool drawelem = false;

            Point<3> newp;
            while(!stepper.GetNextData(newp,dummyt,h) && elnr != -1)
              {
                elnr = GetVolElementHinted(mesh, newp, lami, elnr);
                if(elnr != -1)
                  {
                    drawelem = func(elnr, lami, v);
                    if(dir == -1) v *= -1.;
                    stepper.FeedNextF(v);
                  }
              }

            if (elnr == -1)
              {
                //cout << "direction " <<dir << " reached the wall." << endl;
                break;
              }

            points.Append(newp);
            vals.Append(v.Length());
            drawelems.Append(drawelem);

            if(maxpoints > 0 && points.Size() >= maxpoints)
              {
                break;
              }

            //cout << "length " << length << " h " << h << " vals.Last() " << vals.Last()  << " maxlength " << maxlength << endl;
          }
        dirstart.Append(points.Size());
      }
  }
  
}
