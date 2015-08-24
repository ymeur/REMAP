#include "intersection_ym.hpp"
#include "elt.h"
#include "gpc.h"
#include "grid.h"
#include "triple.h"
#include "polyg.h"
#include <vector>
#include "stdlib.h"

#define epsilon 1e-3  // epsilon distance ratio over side lenght for approximate small circle by great circle
#define fusion_vertex 1e-15 

using namespace std;

void intersect_ym(Elt *a, Elt *b)
{
 
  vector<Coord> srcPolygon ;
  createGreatCirclePolygon(*b, srcGrid.pole, srcPolygon) ;
  vector<Coord> dstPolygon ;
  createGreatCirclePolygon(*a, tgtGrid.pole, dstPolygon) ;


  int na=dstPolygon.size() ;
  Coord *a_gno   = new Coord[na];
  int nb=srcPolygon.size() ;
  Coord *b_gno   = new Coord[nb];

  Coord OC=barycentre(a->vertex,a->n) ;
  Coord Oz=OC ;
  Coord Ox=crossprod(Coord(0,0,1),Oz) ;
  Ox=Ox*(1./norm(Ox)) ;
  Coord Oy=crossprod(Oz,Ox) ;
  double cos_alpha;
  
  for(int n=0; n<na;n++)
  {
    cos_alpha=scalarprod(OC,dstPolygon[n]) ;
    a_gno[n].x=scalarprod(dstPolygon[n],Ox)/cos_alpha ;
    a_gno[n].y=scalarprod(dstPolygon[n],Oy)/cos_alpha ;
    a_gno[n].z=scalarprod(dstPolygon[n],Oz)/cos_alpha ;
  }
  
  for(int n=0; n<nb;n++)
  {
    cos_alpha=scalarprod(OC,srcPolygon[n]) ;
    b_gno[n].x=scalarprod(srcPolygon[n],Ox)/cos_alpha ;
    b_gno[n].y=scalarprod(srcPolygon[n],Oy)/cos_alpha ;
    b_gno[n].z=scalarprod(srcPolygon[n],Oz)/cos_alpha ;
//    cout <<"polygonPoints.InsertPoint("<<n+na<<", "<<b_gno[n].x<<", "<<b_gno[n].y<<", "<<b_gno[n].z<<")"<<endl ; 
  }
//  cout<<"**********************************************"<<endl ;
  
  gpc_polygon src ;
  src.num_contours=1 ;
  src.hole = 0 ;
  src.contour = (gpc_vertex_list*) malloc(sizeof(gpc_vertex_list)) ;
  src.contour->num_vertices = na ;
  src.contour->vertex = (gpc_vertex*) malloc(sizeof(gpc_vertex)*na) ;
  for(int n=0; n<na;n++) 
  {
    src.contour[0].vertex[n].x=a_gno[n].x ;
    src.contour[0].vertex[n].y=a_gno[n].y ;
  }
  
  gpc_polygon trg ;
  trg.num_contours=1 ;
  trg.hole = 0 ;
  trg.contour = (gpc_vertex_list*) malloc(sizeof(gpc_vertex_list)) ;
  trg.contour->num_vertices = nb ;
  trg.contour->vertex = (gpc_vertex*) malloc(sizeof(gpc_vertex)*nb) ;
  for(int n=0; n<nb;n++) 
  {
    trg.contour[0].vertex[n].x=b_gno[n].x ;
    trg.contour[0].vertex[n].y=b_gno[n].y ;
  }

  gpc_polygon intersection ;
  gpc_polygon_clip( GPC_INT,&src,&trg,&intersection) ;
  
//  cout<<"**********************************************"<<endl ;
//  cout<<"Intersection"<<endl ;
//  cout<<intersection.num_contours<<endl ;
  
//  for(int nc=0; nc<intersection.num_contours;nc++)
//   for(int n=0; n < intersection.contour[nc].num_vertices; n++) 
//   {
//    cout<<intersection.contour[nc].vertex[n].x << ", "<< intersection.contour[nc].vertex[n].y <<endl ;
//   }
//  cout<<"**********************************************"<<endl ;

  if (intersection.num_contours==1)
  {
    Coord* intersectPolygon=new Coord[intersection.contour[0].num_vertices] ; 
    for(int n=0; n < intersection.contour[0].num_vertices; n++) 
    {
      intersectPolygon[n]=Ox*intersection.contour[0].vertex[n].x+Oy*intersection.contour[0].vertex[n].y+Oz ;
      intersectPolygon[n]=intersectPolygon[n]*(1./norm(intersectPolygon[n])) ;
//      cout <<"polygonPoints.InsertPoint("<<n<<", "<<intersectPolygon[n].x<<", "<<intersectPolygon[n].y<<", "<<intersectPolygon[n].z<<")"<<endl ;
    }

// remove redondants vertex
    int nv=0 ;
    for(int n=0; n < intersection.contour[0].num_vertices; n++) 
    {
		  if (norm(intersectPolygon[n]-intersectPolygon[(n+1)%intersection.contour[0].num_vertices])>fusion_vertex)
		  {
        intersectPolygon[nv]=intersectPolygon[n] ;
        nv++ ;
      }
    }
//    cout<<"**********************************************"<<endl ;
//    cout << "number of vertex "<<nv<<endl ;
    
//    if (nv<=2) cout<<"Interesection polygon area "<<0<<endl ; 
//    else cout<<"Interesection polygon area : "<<polygonarea(intersectPolygon,intersection.contour[0].num_vertices)<<endl ; 
//    cout<<"**********************************************"<<endl ;
    
    if (nv>2) 
    {
   		Polyg *is = new Polyg;
  		is->x = exact_barycentre(intersectPolygon,nv);
  		is->area = polygonarea(intersectPolygon,nv) ;
  		is->id = b->id; /* intersection holds id of corresponding source element (see Elt class definition for details about id) */
  		is->src_id = b->src_id;
  		is->n = nv;
  		(a->is).push_back(is);
  		(b->is).push_back(is);
    }
    delete[] intersectPolygon ;
  }
  gpc_free_polygon(&src);
  gpc_free_polygon(&trg);  
  gpc_free_polygon(&intersection);  
  delete[] a_gno ;
  delete[] b_gno ;
}    

void createGreatCirclePolygon(const Elt& element, const Coord& pole, vector<Coord>& coordinates)
{
  int nv = element.n; 

  double z,r ;
  int north ;
  int iterations ;
  
  Coord xa,xb,xi,xc ;
  Coord x1,x2,x ;

  for(int i=0;i < nv ;i++)
  {
    north = (scalarprod(element.edge[i], pole) < 0) ? -1 : 1;
    z=north*element.d[i] ;

    if (z != 0.0)
    {

      xa=element.vertex[i] ;
      xb=element.vertex[(i+1)%nv] ;
      iterations=0 ;

// compare max distance (at mid-point) between small circle and great circle
// if greater the epsilon refine the small circle by dividing it recursively.

      do
      {
        xc = pole * z ;
        r=sqrt(1-z*z) ;
        xi=(xa+xb)*0.5 ;
        x1=xc+(xi-xc)*(r/norm(xi-xc)) ;
        x2= xi*(1./norm(xi)) ;
        ++iterations;
        xb=x1 ;
      } while(norm(x1-x2)/norm(xa-xb)>epsilon) ;
      
      iterations = 1 << (iterations-1) ;
      
// small circle divided in "iterations" great circle arc      
      Coord delta=(element.vertex[(i+1)%nv]-element.vertex[i])*(1./iterations);
      x=xa ;
      for(int j=0; j<iterations ; j++)
      {
        //xc+(x-xc)*r/norm(x-xc)
        coordinates.push_back(xc+(x-xc)*(r/norm(x-xc))) ;
        x=x+delta ;
      }
    }
    else coordinates.push_back(element.vertex[i]) ;
  }
}

