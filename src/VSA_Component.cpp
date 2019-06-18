///////////////////////////////////////////////////////////////////////////
// Author: Guillaume Lavoué
// Year: 2008
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
//
//According to: Variational Shape Approximation
//		David Cohen-Steiner, Pierre Alliez and Mathieu Desbrun, SIGGRAPH '2004.
///////////////////////////////////////////////////////////////////////////
#include <mepp_config.h>
#ifdef BUILD_component_VSA

#ifdef _MSC_VER
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#endif

#include "VSA_Component.h"
#include <CGAL/Timer.h>
#include<set>
#include<map>
#include<list>
#include<vector>
#include<stack>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/basic.h>
#include <CGAL/OpenNL/linear_solver.h>
#include <iterator>
using namespace CGAL;

const static  double m_pi_val = 3.14159265359;
const static  double minWeightLap = 0.000001;
const static  double maxWeightLap = 500;

typedef CGAL::Simple_cartesian<double>					AABB_Kernel;
typedef AABB_Kernel::Point_3                            Point;
//typedef CGAL::Vector_3<AABB_Kernel>                     Vector;
typedef CGAL::AABB_polyhedron_triangle_primitive<AABB_Kernel,Polyhedron> AABB_Primitive;
typedef CGAL::AABB_traits<AABB_Kernel, AABB_Primitive>								AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits>												AABB_Tree;
typedef AABB_Tree::Object_and_primitive_id Object_and_primitive_id;
typedef AABB_Tree::Point_and_primitive_id Point_and_primitive_id;

typedef AABB_Kernel::Segment_3 Segment;
typedef AABB_Kernel::Ray_3 Ray;

typedef OpenNL::DefaultLinearSolverTraits <double>  SparseLA;
typedef SparseLA::Vector Vector_OpenNL;
typedef SparseLA::Matrix Matrix;

 typedef struct
  {
	Facet_iterator Facet;

	double DistanceLLoyd;//distance for the LLoyd algorithm

	double PossibleCluster;//cluster for the LLoyd algorithm

  }
  FacetToIntegrate;

   typedef struct
  {
	Vertex_iterator Vertex;
	std::vector<int> TabAdjAnch;


  }
  AnchorVertex;

  double AreaFacetTriangleSeg(Facet_iterator &f)
	{
		Halfedge_around_facet_circulator pHalfedge = f->facet_begin();
		Point3d P = pHalfedge->vertex()->point();
		Point3d Q = pHalfedge->next()->vertex()->point();
		Point3d R = pHalfedge->next()->next()->vertex()->point();

		Vector PQ=Q-P;
                //Vector PR=R-P; // MT
		Vector QR=R-Q;


		Vector normal	=	CGAL::cross_product(PQ,QR);
		double area=0.5*sqrt(normal*normal);

		return area;

	}


	void VSA_Component::Get_Clicked_Vertices(PolyhedronPtr pMesh, double x, double y, int tolerance)
{
	GLdouble *model ;  GLdouble *proj ;  GLint *view;

	view=new GLint[4096];
	model=new double[4096];
	proj=new double[4096];

	glGetIntegerv (GL_VIEWPORT, view);
	glGetDoublev (GL_MODELVIEW_MATRIX, model);
	glGetDoublev (GL_PROJECTION_MATRIX, proj);

	y=view[3]-y;

	GLdouble wx ; GLdouble wy ; GLdouble wz;

	int vertexID=0;
	for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
		gluProject (pVertex->point().x(),pVertex->point().y(),pVertex->point().z(), model, proj, view, &wx, &wy, &wz);  // on simule la projection du sommet dans l'espace window
		if (wz>0. && wz<1)
		if (x>floor(wx)-tolerance && x<floor(wx)+tolerance)
		if (y>floor(wy)-tolerance && y<floor(wy)+tolerance)  // on fixe une petite tol�rance (de 2*5 pixels par exemple) pour la s�lection
		{
		    // Read depth buffer at 2D coordinates obtained from above
            GLfloat bufDepth = 0.0;
            glReadPixels(   static_cast<GLint>( wx ), static_cast<GLint>( wy ),     // Cast 2D coordinates to GLint
                            1, 1,                                                               // Reading one pixel
                            GL_DEPTH_COMPONENT, GL_FLOAT,
                            &bufDepth);

            // Compare depth from buffer to 2D coordinate "depth"
            GLdouble EPSILON = 0.0001;  // Define your own epsilon
            if (fabs(bufDepth - wz) < EPSILON)
            {
                cout<<"3D point is not occluded ";
                pVertex->color(1., 0., 0.);
                cout<<"Vertex ID: "<<vertexID<<endl;
                pMesh->VertexIDSelected = vertexID;
            }
			//statusString.Printf(_T("Vertex: %u  -  (%f, %f, %f)"), vertexID, pVertex->point().x(), pVertex->point().y(), pVertex->point().z());
			//pFrame->set_status_message(statusString);
		}
		vertexID++;
	}

	delete[]view;
	delete[]model;
	delete[]proj;
}

bool sphere_clip_vectorVSA(Point3d &O, double r,const Point3d &P, Vector &V)
{

        Vector W = P - O ;
        double a = (V*V);
        double b = 2.0 * V * W ;
        double c = (W*W) - r*r ;
        double delta = b*b - 4*a*c ;
        if (delta < 0) {
            // Should not happen, but happens sometimes (numerical precision)
            return true ;
        }
        double t = (- b + ::sqrt(delta)) / (2.0 * a) ;
        if (t < 0.0) {
            // Should not happen, but happens sometimes (numerical precision)
            return true ;
        }
        if (t >= 1.0) {
            // Inside the sphere
            return false ;
        }

        V=V*t;

        return true ;
}


void sphere_search(PolyhedronPtr pMesh, Vertex *pVertex, double radius)
{
    Point3d O = pVertex->point() ;

	for (Vertex_iterator pV = pMesh->vertices_begin(); pV!= pMesh->vertices_end(); pV++)
	{
        Point3d P = pV->point() ;

        Point3d p1 = pVertex->point() ;
        Point3d p2 = pV->point();
        Vector V = (p2-p1);

        bool isect = sphere_clip_vectorVSA(O, radius, P, V) ;
        if (!isect) {
            pV->tagged_vertex = 1;

        }
	}
}

void geodes_sphere_search(Vertex *pVertex, double radius)
{

       std::set<Vertex*> vertices ;
        Point3d O = pVertex->point() ;
        std::stack<Vertex*> S ;
        S.push(pVertex) ;
        vertices.insert(pVertex) ;
		int iter=0;
        while(!S.empty())
		{
			Vertex* v = S.top() ;
            S.pop() ;
            Point3d P = v->point() ;
            Halfedge_around_vertex_circulator h = v->vertex_begin();
			Halfedge_around_vertex_circulator pHalfedgeStart = h;
			CGAL_For_all(h,pHalfedgeStart)
			{
                if (!h->is_border_edge())
                {
                    Point3d p1 = h->vertex()->point();
                    Point3d p2 = h->opposite()->vertex()->point();
                    Vector V = (p2-p1);
    //                if (v==pVertex || V * (P - O) > 0.0)
    //				{
                                                //pVertex->color(0., 1., 0.);

                        //double len_old = std::sqrt(V*V);
                        bool isect = sphere_clip_vectorVSA(O, radius, P, V) ;
                        if (!isect) {
                            Vertex_iterator w=h->opposite()->vertex();
    //						h->opposite()->vertex()->color(0., 1., 0.);
                            h->opposite()->vertex()->tagged_vertex = 1;
                            if (vertices.find(&(*w)) == vertices.end()) {
                                vertices.insert(&(*w)) ;
                                S.push(&(*w)) ;
                            }
                        }
    //				}
                }
			}
			iter++;
		}
		double area=m_pi_val*radius*radius;
}



void VSA_Component::draw_Sphere(PolyhedronPtr pMesh, float r, int sphereDisk)
{

    cout<<"ok1"<<endl;
    Timer timer;
	timer.start();

    Vertex_iterator pVertex = NULL;
	for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
        pVertex->tagged_vertex = 0;
        pVertex->convex_component = 0;
        pVertex->visited = 0;
        pVertex->color(0.5, 0.5, 0.5);

	}

    int k=0;
    pVertex = NULL;
	for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
        if(k==pMesh->VertexIDSelected)
        {
            if(sphereDisk == 0)
            sphere_search(pMesh, (&(*pVertex)), r);
            else
            geodes_sphere_search((&(*pVertex)), r);
        }

        k+=1;
    }

//    pVertex = NULL;
//	for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
//	{
//        if(pVertex->tagged_vertex ==1)
//        {
//            pVertex->color(0., 1., 0.);
//            cout<<"tagged vertex detected"<<endl;
//        }
//	}

	//detect border of intersection regions
    pVertex = NULL;
	for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
        if(pVertex->tagged_vertex ==1)
        {
            pVertex->color(0., 0., 1.);
            Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
            Halfedge_around_vertex_circulator pHalfedgeStart = h;
            CGAL_For_all(h,pHalfedgeStart)
            {
				if(h->opposite()->vertex()->tagged_vertex == 0)
                    pVertex->tagged_vertex = 2;
            }
        }
	}

    std::stack<Vertex*> verticesInBorder ;
    pVertex = NULL;
	for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
        if(pVertex->tagged_vertex ==2)
        {
            verticesInBorder.push((&(*pVertex)));
            pVertex->color(0., 1., 0.);
            //cout<<"tagged vertex detected"<<endl;
        }
	}

    int component=0;
    while(!verticesInBorder.empty())
    {
        if(verticesInBorder.top()->convex_component == 0)
        {
            component += 1;
            verticesInBorder.top()->convex_component = component;
        }
//        verticesInBorder.top()->visited = 1;
        Halfedge_around_vertex_circulator h = verticesInBorder.top()->vertex_begin();
        Halfedge_around_vertex_circulator pHalfedgeStart = h;
                verticesInBorder.pop() ;

        CGAL_For_all(h,pHalfedgeStart)
        {
			if( (h->opposite()->vertex()->tagged_vertex == 2) )
            {
                if((h->opposite()->vertex()->convex_component) == 0)
                {
                    h->opposite()->vertex()->convex_component = component;
                    verticesInBorder.push((&(*h->opposite()->vertex())));
                }
//                else
//                    h->opposite()->vertex()->convex_component = h->opposite()->vertex()->convex_component;
//                h->opposite()->vertex()->visited = 1;
            }
        }

    }




    pVertex = NULL;
	for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
        if(pVertex->convex_component != 0)
        {
            double R=(double)(pVertex->convex_component)/(double)component*255.;
			int indiceLut=floor(R);
            pVertex->color(LUT_Seg[3*(indiceLut)],LUT_Seg[3*(indiceLut)+1],LUT_Seg[3*(indiceLut)+2]);
        }
        if(k==pMesh->VertexIDSelected)
        {
            pVertex->color(1., 0., 0.);
        }
        k+=1;
//        if(pVertex->convex_component ==1)
//        {
//            pVertex->color(1., 0., 0.);
//        }
	}
    cout<<"connex components number" <<component<<"  comput. time = "<<timer.time()<<"seconds"<<endl;
}


double VSA_Component::AreaFacetTriangle(Facet_handle &f)
	{
		Halfedge_around_facet_circulator pHalfedge = f->facet_begin();
		Point3d P = pHalfedge->vertex()->point();
		Point3d Q = pHalfedge->next()->vertex()->point();
		Point3d R = pHalfedge->next()->next()->vertex()->point();

		Vector PQ=Q-P;
                //Vector PR=R-P; // MT
		Vector QR=R-Q;


		Vector normal	=	CGAL::cross_product(PQ,QR);
		double area=0.5*sqrt(normal*normal);

		return area;

	}

  void VSA_Component::Init(int NbProxy)
	{
		////Creation of the initial proxies by random seed triangle picking
		m_Table_Proxy.clear();

		m_NbProxy=NbProxy;
		int NbFacet=m_Poly->size_of_facets();
		int offset=NbFacet/NbProxy;
		int i=0;
		for(Facet_iterator	pface	=	m_Poly->facets_begin();
				pface	!= m_Poly->facets_end();
				pface++)
		{

			if(i%offset==0)///this triangle is chosen
			{
				//a proxy is created
				Proxy NewProxy;
				NewProxy.Normal=pface->normal();

				NewProxy.Seed=pface;
				///the proxy is added
				m_Table_Proxy.push_back(NewProxy);


			}

			pface->LabelVSA=-1;
			i++;
		}

		m_Poly->NbFaceLabel=m_NbProxy;

	}

	struct CompFacet
	{
	bool operator()(FacetToIntegrate f1, FacetToIntegrate f2) const
	{

		return(f1.DistanceLLoyd<f2.DistanceLLoyd);

	}
	};


	double VSA_Component::DistorsionError(Facet_iterator f,Proxy p)
	{
		Vector v=f->normal()-p.Normal;
		double nrm=v*v;
		double area=AreaFacetTriangleSeg(f);
		return nrm*area;
	}

std::vector<double> VSA_Component::computeLaplacianCotg2(PolyhedronPtr pMesh)
{

	double sumAlljWeight=0.0 ;
	std::vector<double> ListLaplaceCot;
	for(Vertex_iterator	pVertex	= pMesh->vertices_begin();
		pVertex	!= pMesh->vertices_end();
		pVertex++)
	{

		double L=0.0;
		double sumjweight=0.0;
		Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
		Halfedge_around_vertex_circulator pHalfedgeStart = h;
		CGAL_For_all(h,pHalfedgeStart)
		{
			Vertex_handle j=h->opposite()->vertex();
			double jweight = computeWcot(h);
			if (jweight<minWeightLap)
				jweight = minWeightLap;
			if (jweight>maxWeightLap)
				jweight = maxWeightLap;

			L+= jweight * ( (j->KminCurv * j->KmaxCurv) );
			sumjweight+=jweight;
		}
//        if (sumjweight < minWeightLap)
//        {
        pVertex->LaplaceCot = std::abs(  ( pVertex->KminCurv * pVertex->KmaxCurv) - (L/sumjweight) )  ;
       // cout << "pVertex->KmaxCurv" << pVertex->KmaxCurv << endl;
//        }
//        else
//        {
//         pVertex->LaplaceCot = 0;
//        }

		sumAlljWeight+=pVertex->LaplaceCot;

		ListLaplaceCot.push_back(pVertex->LaplaceCot);
	}

	return ListLaplaceCot;
}

std::vector<double> VSA_Component::computeLaplacianCotg3(PolyhedronPtr pMesh)
{

	double sumAlljWeight=0.0 ;
	std::vector<double> ListLaplaceCot;
	for(Vertex_iterator	pVertex	= pMesh->vertices_begin();
		pVertex	!= pMesh->vertices_end();
		pVertex++)
	{

		double L=0.0;
		double sumjweight=0.0;
		Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
		Halfedge_around_vertex_circulator pHalfedgeStart = h;
		CGAL_For_all(h,pHalfedgeStart)
		{
			Vertex_handle j=h->opposite()->vertex();
			double jweight = computeWcot(h);
			if (jweight<minWeightLap)
				jweight = minWeightLap;
			if (jweight>maxWeightLap)
				jweight = maxWeightLap;

			L+= jweight * ( j->vertGC );
			sumjweight+=jweight;
		}
//        if (sumjweight < minWeightLap)
//        {
        pVertex->LaplaceCot = std::abs(  ( pVertex->vertGC) - (L/sumjweight) )  ;

//        }
//        else
//        {
//         pVertex->LaplaceCot = 0;
//        }

		sumAlljWeight+=pVertex->LaplaceCot;

		ListLaplaceCot.push_back(pVertex->LaplaceCot);
	}

	return ListLaplaceCot;
}

double VSA_Component::computeWcot(Halfedge_handle h)
{
	double res=0.0; double count=0.0;
	Facet_handle pFacet1 = h->facet();
	Facet_handle pFacet2 = h->opposite()->facet();
	CGAL_assertion(pFacet1 != pFacet2);
	if (pFacet1 != NULL)
	{
		res+= computeHalfWcot(h);
		count+=1.0;
	}
	if (pFacet2 != NULL)
	{
		res+= computeHalfWcot(h->opposite());
		count+=1.0;
	}
	return (res/count);
}

double VSA_Component::computeHalfWcot(Halfedge_handle h)
{
	Point3d B=h->next()->vertex()->point();
	Point3d J=h->vertex()->point();
	Point3d I=h->next()->next()->vertex()->point();
	Vector BJ=J-B;
	Vector BI=I-B;
	double CosAngle = (BJ * BI)/(std::sqrt(BJ*BJ) * std::sqrt(BI*BI)) ;
	return ( 1.0 / std::tan(acos_stable(CosAngle)));
}

double VSA_Component::acos_stable(double CosAng)
{
	if (CosAng >=1.0)
		CosAng= 1.0;
	else if(CosAng <= -1.0)
		CosAng= -1.0;

	return (std::acos(CosAng)) ;
}


double VSA_Component::getMaxDim(PolyhedronPtr polyhedron_ptr)
{
	polyhedron_ptr->compute_bounding_box();

	double max=polyhedron_ptr->xmax()-polyhedron_ptr->xmin();
	if(polyhedron_ptr->ymax()-polyhedron_ptr->ymin()>max)
		max = polyhedron_ptr->ymax()-polyhedron_ptr->ymin();
	if(polyhedron_ptr->zmax()-polyhedron_ptr->zmin()>max)
		max = polyhedron_ptr->zmax()-polyhedron_ptr->zmin();

	return max;

}

double VSA_Component::getAngle(Vector he1,Vector he2)

{

    double denum = std::sqrt(he1*he1) * std::sqrt(he2*he2);

    double cosAngleNeighEdges;

    if(denum!=0)

        cosAngleNeighEdges = (he1 * he2) / denum;    //may cause prob in quantization

    else

        return(0.0);

    return acos_stable(cosAngleNeighEdges);

}


void VSA_Component::computeGaussianCurv1Ring(PolyhedronPtr pMesh)

{

    //double maxDimInner = getMaxDimInner(pMesh);

       double MinNrmGauss2Curvature=100000.0;

       double MaxNrmGauss2Curvature=-100000.0;

        for(Vertex_iterator        pVertex        = pMesh->vertices_begin();

                pVertex        != pMesh->vertices_end();

                pVertex++)

            pVertex->vertGC = (double)0.0;

    for( Halfedge_iterator phEdges = pMesh->halfedges_begin() ;

        phEdges != pMesh->halfedges_end();

        phEdges++){

            if (!phEdges->is_border_edge())

           {

            Vector edge1 =  phEdges->opposite()->vertex()->point() - phEdges->vertex()->point();

            Halfedge_handle heNext;

            heNext = phEdges->next();

            Vector edge2 =  heNext->vertex()->point() - heNext->opposite()->vertex()->point();

            double angleValue = getAngle(edge1,edge2);

            phEdges->vertex()->vertGC = phEdges->vertex()->vertGC + angleValue;

            }

        }

    for(Vertex_iterator        pVertex        = pMesh->vertices_begin();

                pVertex        != pMesh->vertices_end();

                pVertex++)

                {

                    pVertex->vertGC = (2.0 * m_pi_val) - pVertex->vertGC;

                    MinNrmGauss2Curvature=CGAL::min<double>(MinNrmGauss2Curvature,pVertex->vertGC);

            MaxNrmGauss2Curvature=CGAL::max<double>(MaxNrmGauss2Curvature,pVertex->vertGC);

                }

}

std::vector<double> VSA_Component::DeviationPolyhedronToPolyhedron(PolyhedronPtr original, PolyhedronPtr lisse)
{

	// constructs AABB tree
    AABB_Tree tree(lisse->facets_begin(),lisse->facets_end());
    tree.accelerate_distance_queries();



	Vertex_iterator	pVertex;
	int numVertex = original->size_of_vertices();
	//Vector * DistanceDeviation = new Vector[numVertex];

	std::vector<double> DistanceDeviation;
	//std::string text= "XPoint1(OriginalMesh) YPoint1(OriginalMesh) ZPoint1(OriginalMesh) XPoint2(SmoothedMesh) YPoint2(SmoothedMesh) ZPoint2(SmoothedMesh) XVectorDeviation YVectorDeviation ZVectorDeviation Distance";

	for (pVertex = original->vertices_begin(); pVertex != original->vertices_end(); pVertex++)
		{


            Point_and_primitive_id pp = tree.closest_point_and_primitive(pVertex->point());
            Point3d Nearest=pp.first;
            Polyhedron::Face_handle f = pp.second;
            Vector V_normal = f->halfedge()->vertex()->normal();
            // query point
			Point query = pVertex->point();
			Vector deviation= query - Nearest;
			double angle= Get_Angle_Weight (deviation,V_normal);
			double aa= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));
			if (angle < 0)
			{
			  aa = -aa;
			}


            DistanceDeviation.push_back(aa);
            pVertex->depthvertex = aa;
	       }

	 WriteTxt2(DistanceDeviation, "/home/isit/dev.txt");

	return DistanceDeviation;
 }

void VSA_Component::WriteTxt2(std::vector<double> dev, std::string file)
{
	std::ofstream f;

	f.open(file.c_str(),std::ios::out|std::ios::trunc);

  if (!f.is_open())

   std::cout << "Impossible d'ouvrir le fichier en ecriture !" << endl;

  else

  {


	  for (int i = 0; i < dev.size(); i++)
	  {

		  f << dev.at(i) << " \n";

	  }
   }

  f.close();

}

PolyhedronPtr VSA_Component::TranslateToZeroTheDepth(PolyhedronPtr pMesh, std::vector<double> dev)
{
  // trouver le minimum du profondeur
  double minimum= dev.at(0);
  int i=1;
  while (i < dev.size())
	  {
           if (minimum > dev.at(i))
           {
                minimum = dev.at(i);
           }
            i++;
	  }

cout << "minimum =" << minimum << endl;
cout << "abs minimum =" << abs(minimum) << endl;
 for(Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)

   {
        pVertex->depthvertex= pVertex->depthvertex + abs(minimum);
   }

  return pMesh;

}

PolyhedronPtr VSA_Component::DepthFacets(PolyhedronPtr pMesh)
{
   double sum =0.0;
   int nbr =0;
   double depth;
    for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	{
        Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
        sum =0.0;
        nbr =0;
		do
		{
			sum= sum + pHalfedge->vertex()->depthvertex;
			nbr =nbr +1;
		}
		while(++pHalfedge != pFacet->facet_begin());
		depth = sum/nbr;
		if ( depth < 0)
		{
		depth= -depth;
		}
		if (isnan(depth))
		{
		depth=0.0;
		}
        pFacet->depthfacet= depth;

	}

	return pMesh;
}



double VSA_Component::getOtsuThreshold(std::vector<double> data)
 {
		// Otsu's threshold algorithm
		// C++ code by Jordan Bevik <Jordan.Bevic@qtiworld.com>
		// ported to ImageJ plugin by G.Landini
		double kStar;  // k = the current threshold; kStar = optimal threshold
		int N1, N;    // N1 = # points with intensity <=k; N = total number of points
		double BCV, BCVmax; // The current Between Class Variance and maximum BCV
		double num, denom;  // temporary bookeeping
		int Sk;  // The total intensity for all histogram points <=k
		int S, L=256; // The total intensity of the image

		// Initialize values:
		S = N = 0;
		int k;
		for (k=0; k<L; k++){
			S += k * data[k];	// Total histogram intensity
			N += data[k];		// Total number of data points
		}

		Sk = 0;
		N1 = data[0]; // The entry for zero intensity
		BCV = 0;
		BCVmax=0;
		kStar = 0.0;

		// Look at each possible threshold value,
		// calculate the between-class variance, and decide if it's a max
		for (k=1; k<L-1; k++) { // No need to check endpoints k = 0 or k = L-1
			Sk += k * data[k];
			N1 += data[k];

			// The float casting here is to avoid compiler warning about loss of precision and
			// will prevent overflow in the case of large saturated images
			denom = (double)( N1) * (N - N1); // Maximum value of denom is (N^2)/4 =  approx. 3E10

			if (denom != 0 ){
				// Float here is to avoid loss of precision when dividing
				num = ( (double)N1 / N ) * S - Sk; 	// Maximum value of num =  255*N = approx 8E7
				BCV = (num * num) / denom;
			}
			else
				BCV = 0;

			if (BCV >= BCVmax){ // Assign the best threshold found so far
				BCVmax = BCV;
				kStar = k;
			}
		}
		// kStar += 1;	// Use QTI convention that intensity -> 1 if intensity >= k
		// (the algorithm was developed for I-> 1 if I <= k.)
		return kStar;
	}

double VSA_Component::MedianofVector (std::vector<double> grades)
{

    sort(grades.begin(), grades.end());
    double median= *(grades.begin()+grades.size()/2); //89

    return median;
}



 std::list<Facet_iterator> VSA_Component::addFacetToOrderedList(Facet_iterator pVertex,  std::list<Facet_iterator> Vertexs )
{
    if (Vertexs.empty())
    {
        Vertexs.push_front(pVertex);
    }
    else
    {

          std::list<Facet_iterator> ::iterator it = Vertexs.begin();
          Facet_iterator v;
          while (it != Vertexs.end())
          {
            v= *it;
           // cout << "v->gradient_depth " << v->gradient_depth << endl;
            if (v->depthfacet < pVertex->depthfacet) ++it;
             else
             break;
          }
          Vertexs.insert(it,pVertex);

    }

    return Vertexs;

}

PolyhedronPtr VSA_Component::Morpho_Fill_holes_Using_Sphere(PolyhedronPtr pMesh, int indR, double R, double R2)
{
     std::vector<Facet_handle> facets= m_Table_Proxy[indR].facets;
    int  ind= facets[0]->LabelVSA ;
    Facet_iterator pFacet, pFacet2;

    cout << "entrer" << endl;
    for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	    {


            if (pFacet->LabelVSA == ind)
            {
                pFacet->morphofilled=1;
                pFacet->visitedMF =0;
            }

            else
            {
             pFacet->morphofilled=0;
             pFacet->visitedMF =0;
            }

		}

		// trouver tous les facettes qui appartiennent la différence de deux spheres ciconscrites a la region de centre le baryucentre de la région
		// et de rayon R= sqrt de la surface de la région et R2= 2*R

        std::list<Facet_handle> facets_voisines_2;
        Facet_handle dist_facet, faceRef;
        // the barycenter of a region
        Point barycenterR, barycenterF ;
        double dist;
        compute_region_center(pMesh,facets, barycenterR);

         for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
                {
                           pMesh->compute_facet_center(pFacet,barycenterF);
                           Vector deviation= barycenterR - barycenterF;
                           dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));
                           if (dist> R && dist< R2)
                           {
                                    facets_voisines_2.push_back(pFacet);
                                    pFacet->morphofilled= 2;
                           }

                    }

        Facet_handle faceV;

        while(!facets_voisines_2.empty())
        {

             faceV= facets_voisines_2.front();
             facets_voisines_2.pop_front();

             Halfedge_around_facet_circulator pHalfedge = faceV->facet_begin();

            do
            {
                    pMesh->compute_facet_center(pHalfedge->opposite()->facet(),barycenterF);
                    Vector deviation= barycenterR - barycenterF;
                    dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));
                   if ((pHalfedge->opposite()->facet()->morphofilled != 1) && (pHalfedge->opposite()->facet()->visitedMF == 0) && (dist< R2) )
                   {
                    pHalfedge->opposite()->facet()->visitedMF = 1;

                    if (pHalfedge->opposite()->facet()->morphofilled == 0) {
                        pHalfedge->opposite()->facet()->morphofilled= 2;
                        }

                        facets_voisines_2.push_back(pHalfedge->opposite()->facet());

                   }

            }
            while(++pHalfedge != faceV->facet_begin());

        }

 return pMesh;

}

void VSA_Component::compute_region_center(PolyhedronPtr pMesh, std::vector<Facet_handle> facets, Point& center)
	{
		int degree = 0;
        Vector vec(0.0,0.0,0.0);
        Point barycenter;
			for (int i=0; i<facets.size(); i++)
			{
                pMesh->compute_facet_center(facets.at(i),barycenter);
				vec = vec + (barycenter-CGAL::ORIGIN);
				degree++;
			}
			center = CGAL::ORIGIN + (vec/(double)degree);


	}

PolyhedronPtr VSA_Component::Morpho_Fill_holes(PolyhedronPtr pMesh, int indR)
{
    std::vector<Facet_handle> facets= m_Table_Proxy[indR].facets;
    int  ind= facets[0]->LabelVSA ;
    Facet_iterator pFacet, pFacet2;

    cout << "entrer" << endl;
    for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	    {


            if (pFacet->LabelVSA == ind)
            {
                pFacet->morphofilled=1;
                pFacet->visitedMF =0;
            }

            else
            {
             pFacet->morphofilled=0;
             pFacet->visitedMF =0;
            }

		}



		// find a facet the most far away from our region
        double dist, max_dist=0;
        Point barycenter1, barycenter2;
        Facet_handle dist_facet, faceRef;
        //faceRef= m_Table_Proxy[indR].Seed;
        faceRef= facets[0];

        pMesh->compute_facet_center(faceRef,barycenter1);

         for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
                {
                    if (pFacet != faceRef && pFacet->morphofilled==0 )
                    {
                           pMesh->compute_facet_center(pFacet,barycenter2);
                           Vector deviation= barycenter1 - barycenter2;
                           dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));
                           if (max_dist < dist)
                           {
                            max_dist= dist;
                            dist_facet= pFacet;

                           }

                    }
                }

        cout << "max dist ==" << max_dist << endl;
        cout << "dist_facet->filled" << dist_facet->morphofilled << endl;

        std::list<Facet_handle> facets_voisines_2;

        facets_voisines_2.push_back(dist_facet);
        dist_facet->morphofilled= 2;
        Facet_handle faceV;

        while(!facets_voisines_2.empty())
        {

             faceV= facets_voisines_2.front();
             facets_voisines_2.pop_front();

             Halfedge_around_facet_circulator pHalfedge = faceV->facet_begin();

            do
            {

                   if ((pHalfedge->opposite()->facet()->morphofilled != 1) && (pHalfedge->opposite()->facet()->visitedMF == 0) )
                   {
                    pHalfedge->opposite()->facet()->visitedMF = 1;

                    if (pHalfedge->opposite()->facet()->morphofilled == 0) {
                        pHalfedge->opposite()->facet()->morphofilled= 2;
                        }

                        facets_voisines_2.push_back(pHalfedge->opposite()->facet());

                   }

            }
            while(++pHalfedge != faceV->facet_begin());

        }


        //
		return pMesh;

}


void VSA_Component::Erosion(PolyhedronPtr pMesh, int indR, double R )
{


    std::vector<Facet_handle> facets= m_Table_Proxy[indR].facets;
    int  ind= facets[0]->LabelVSA;

cout << "entrer" << endl;
    for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	    {

            pFacet->visited=0;
            if (pFacet->LabelVSA == ind)
            {

                pFacet->eroded=1;
            }

            else
            {
             pFacet->eroded=0;

            }

		}
		cout << "entrer" << endl;

// trouver les facettes qui entourent la région
     Facet_handle face, faceV, pfacet2;

    std::vector<Facet_handle> facets_voisines;


      for (int i=0; i<facets.size(); i++)
        {
                face = facets.at(i);

                Halfedge_around_facet_circulator pHalfedge = face->facet_begin();
                Halfedge_around_facet_circulator end = pHalfedge;

				CGAL_For_all(pHalfedge,end)
				{
					Facet_handle pfacet2 = pHalfedge->opposite()->facet();
					if ((pfacet2 != NULL) && (pfacet2->eroded == 0))
                    {
					                 facets_voisines.push_back(pfacet2);
					}
				}


           }


cout << "size facets_voisines" << facets_voisines.size() << endl;
 std::list<Facet_handle> facets_voisines_2;
Point barycenter1, barycenter2;
double dist;

  for (int i=0; i<facets_voisines.size(); i++)
     {

        face = facets_voisines[i];
        facets_voisines_2.push_back(face);
        pMesh->compute_facet_center(face,barycenter1);
        face->eroded= 0;
        while(!facets_voisines_2.empty())
        {

             faceV= facets_voisines_2.front();
             facets_voisines_2.pop_front();

             Halfedge_around_facet_circulator pHalfedge = faceV->facet_begin();

            do
            {
                  // barycentric subdivision
                   pMesh->compute_facet_center(pHalfedge->opposite()->facet(),barycenter2);
                   Vector deviation= barycenter1 - barycenter2;
                   dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));


                   if (dist < R && pHalfedge->opposite()->facet()->visited == 0 )
                   {
                    pHalfedge->opposite()->facet()->visited = 1;

                    if (pHalfedge->opposite()->facet()->eroded == 1) {
                        pHalfedge->opposite()->facet()->eroded= 0;
                        }

                        facets_voisines_2.push_back(pHalfedge->opposite()->facet());

                   }

            }
            while(++pHalfedge != faceV->facet_begin());

        }


        face = facets_voisines[i];
        facets_voisines_2.push_back(face);
        pMesh->compute_facet_center(face,barycenter1);

        while(!facets_voisines_2.empty())
        {

             faceV= facets_voisines_2.front();
             facets_voisines_2.pop_front();

             Halfedge_around_facet_circulator pHalfedge = faceV->facet_begin();

            do
            {
                  // barycentric subdivision
                   pMesh->compute_facet_center(pHalfedge->opposite()->facet(),barycenter2);
                   Vector deviation= barycenter1 - barycenter2;
                   dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));
                   //cout << "dist" << dist << endl;

                   if (dist < R && pHalfedge->opposite()->facet()->visited == 1)
                   {
                   pHalfedge->opposite()->facet()->visited = 0;
                   facets_voisines_2.push_back(pHalfedge->opposite()->facet());

                   }


            }
            while(++pHalfedge != faceV->facet_begin());

        }


     }




  }

PolyhedronPtr VSA_Component::dilatation(PolyhedronPtr pMesh, int indR, double R )
{


    std::vector<Facet_handle> facets= m_Table_Proxy[indR].facets;
    int  ind= facets[0]->LabelVSA;

   //cout << "entrer" << endl;
    for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	    {


            if (pFacet->LabelVSA == ind)
            {
                pFacet->visited=0;
                pFacet->morpho_dilated=1;
            }

            else
            {
             pFacet->morpho_dilated=0;
             pFacet->visited=0;
            }

		}
		//cout << "entrer" << endl;

// trouver les facettes qui entourent la région
     Facet_handle face, faceV, pfacet2;

    std::vector<Facet_handle> facets_voisines;


      for (int i=0; i<facets.size(); i++)
        {
                face = facets.at(i);

                Halfedge_around_facet_circulator pHalfedge = face->facet_begin();
                Halfedge_around_facet_circulator end = pHalfedge;

				CGAL_For_all(pHalfedge,end)
				{
					Facet_handle pfacet2 = pHalfedge->opposite()->facet();
					if ((pfacet2 != NULL) && (pfacet2->morpho_dilated == 0) && (face->visited==0))
                    {
					                 facets_voisines.push_back(face);
					                 face->visited= 1;
					}
				}


           }


       for (int i=0; i<facets.size(); i++)
            {
              facets.at(i)->visited=0;
            }

//cout << "size facets_voisines" << facets_voisines.size() << endl;
 std::list<Facet_handle> facets_voisines_2;
Point barycenter1, barycenter2;
double dist;

  for (int i=0; i<facets_voisines.size(); i++)
     {
        //std::cout << "on reprend avec une nouvelle facette: " << i << std::endl;
        face = facets_voisines[i];
        facets_voisines_2.push_back(face);
        pMesh->compute_facet_center(face,barycenter1);
        face->morpho_dilated= 1;
        while(!facets_voisines_2.empty())
        {

             faceV= facets_voisines_2.front();
             facets_voisines_2.pop_front();

             Halfedge_around_facet_circulator pHalfedge = faceV->facet_begin();

            do
            {
                  // barycentric subdivision
                   pMesh->compute_facet_center(pHalfedge->opposite()->facet(),barycenter2);
                   Vector deviation= barycenter1 - barycenter2;
                   dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));
                   //cout << "dist" << dist << endl;

                   if (dist < R && pHalfedge->opposite()->facet()->visited == 0 )
                   {
                    pHalfedge->opposite()->facet()->visited = 1;

                    if (pHalfedge->opposite()->facet()->morpho_dilated == 0) {
                        pHalfedge->opposite()->facet()->morpho_dilated= 1;
                        }

                        facets_voisines_2.push_back(pHalfedge->opposite()->facet());
                      //  std::cout << " on supprime une facette" << std::endl;
                   }

            }
            while(++pHalfedge != faceV->facet_begin());

        }

       // std::cout << "on nettoie" << std::endl;
        face = facets_voisines[i];
        facets_voisines_2.push_back(face);
        pMesh->compute_facet_center(face,barycenter1);
        while(!facets_voisines_2.empty())
        {

             faceV= facets_voisines_2.front();
             facets_voisines_2.pop_front();

             Halfedge_around_facet_circulator pHalfedge = faceV->facet_begin();

            do
            {
                  // barycentric subdivision
                   pMesh->compute_facet_center(pHalfedge->opposite()->facet(),barycenter2);
                   Vector deviation= barycenter1 - barycenter2;
                   dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));
                   //cout << "dist" << dist << endl;

                   if (dist < R && pHalfedge->opposite()->facet()->visited == 1)
                   {
                   pHalfedge->opposite()->facet()->visited = 0;
                   facets_voisines_2.push_back(pHalfedge->opposite()->facet());

                   }


            }
            while(++pHalfedge != faceV->facet_begin());

        }
        //std::cout << "on passe à un triangle de la mer suivant" << std::endl;

     }
       //f std::cout << "on a fini" << std::endl;



return pMesh;

  }

PolyhedronPtr VSA_Component::Erosion_after_dilation(PolyhedronPtr pMesh, int indR, double R )
{
    std::vector<Facet_handle> facets;
    for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	    {

                pFacet->morpho_closed = pFacet->morpho_dilated;


                if (pFacet->morpho_closed != 0)
                facets.push_back(pFacet);

                pFacet->visited= 0;

		}

     // trouver les facettes qui entourent la région
    // trouver les facettes qui entourent la région
     Facet_handle face, faceV, pfacet2;

    std::vector<Facet_handle> facets_voisines;


      for (int i=0; i<facets.size(); i++)
        {
                face = facets.at(i);

                Halfedge_around_facet_circulator pHalfedge = face->facet_begin();
                Halfedge_around_facet_circulator end = pHalfedge;

				CGAL_For_all(pHalfedge,end)
				{
					Facet_handle pfacet2 = pHalfedge->opposite()->facet();
					if ((pfacet2 != NULL) && (pfacet2->morpho_closed == 0) && (pfacet2->visited == 0) )
                    {
					   facets_voisines.push_back(pfacet2);
					   pfacet2->visited = 1;
					}
				}

           }

       for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	    {

               pFacet->visited= 0;
        }

 //cout << "size facets_voisines" << facets_voisines.size() << endl;
 std::list<Facet_handle> facets_voisines_2;
Point barycenter1, barycenter2;
double dist;

  for (int i=0; i<facets_voisines.size(); i++)
     {

        face = facets_voisines[i];
        facets_voisines_2.push_back(face);
        pMesh->compute_facet_center(face,barycenter1);
        face->morpho_closed= 0;
        while(!facets_voisines_2.empty())
        {

             faceV= facets_voisines_2.front();
             facets_voisines_2.pop_front();

             Halfedge_around_facet_circulator pHalfedge = faceV->facet_begin();

            do
            {
                  // barycentric subdivision
                   pMesh->compute_facet_center(pHalfedge->opposite()->facet(),barycenter2);
                   Vector deviation= barycenter1 - barycenter2;
                   dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));


                   if (dist < R && pHalfedge->opposite()->facet()->visited == 0 )
                   {
                    pHalfedge->opposite()->facet()->visited = 1;

                    if (pHalfedge->opposite()->facet()->morpho_closed == 1) {
                        pHalfedge->opposite()->facet()->morpho_closed= 0;
                        }

                        facets_voisines_2.push_back(pHalfedge->opposite()->facet());

                   }

            }
            while(++pHalfedge != faceV->facet_begin());

        }

        facets_voisines_2.push_back(face);

        while(!facets_voisines_2.empty())
        {

             faceV= facets_voisines_2.front();
             facets_voisines_2.pop_front();

             Halfedge_around_facet_circulator pHalfedge = faceV->facet_begin();

            do
            {
                  // barycentric subdivision
//                   pMesh->compute_facet_center(pHalfedge->opposite()->facet(),barycenter2);
//                   Vector deviation= barycenter1 - barycenter2;
//                   dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));
                   //cout << "dist" << dist << endl;

                   if (pHalfedge->opposite()->facet()->visited == 1)
                   {
                   pHalfedge->opposite()->facet()->visited = 0;
                   facets_voisines_2.push_back(pHalfedge->opposite()->facet());

                   }


            }
            while(++pHalfedge != faceV->facet_begin());

        }


     }

		return pMesh;
}

PolyhedronPtr VSA_Component::closing(PolyhedronPtr pMesh, int indR, double R )
{
    //dilatation suivie d'une érosion
    pMesh= dilatation(pMesh,indR,R );
    //cout << "end dilatation " << endl;
    pMesh= Erosion_after_dilation(pMesh,indR,R);

    return pMesh;
}

void VSA_Component::ReorganizeListRegions(PolyhedronPtr pMesh)
{
    int j=1;
    double area;
    for(int i=0;i<m_Table_Proxy.size();i++)//For each proxy we select triangles that could grow
	{
         m_Table_Proxy[i].Label= j;
         area=0.0;
         for (int k=0;k<m_Table_Proxy[i].facets.size();k++)
         {
            m_Table_Proxy[i].facets[k]->LabelVSA= j;
            area = area + AreaFacetTriangle(m_Table_Proxy[i].facets[k]);
         }
         m_Table_Proxy[i].Area= area;
         m_Table_Proxy[i].Normal= m_Table_Proxy[i].Seed->normal();
         j++;
	}
	//cout << "j="<< j << "nber of regions =" << pMesh->NbFaceLabel << endl;
}




bool VSA_Component::distFtoList(PolyhedronPtr pMesh, Facet_handle face, std::vector<Facet_handle> facets, double R)
{
   bool exist = false;
   int i=0;
   Point barycenter1, barycenter2;
   double dist;

   pMesh->compute_facet_center(face,barycenter1);
   while (!exist && i < facets.size())
   {
     pMesh->compute_facet_center(facets.at(i),barycenter2);
     Vector deviation= barycenter1 - barycenter2;
     dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));
    // cout << "dist" << dist << endl;
      if (dist > R)
      {
            exist= true;
           // cout << "yeeees" << endl;

      }


      i=i+1;
   }
   return exist;
}



PolyhedronPtr VSA_Component::Burst_Wind_Segmentation_ON_Facets_2(PolyhedronPtr pMesh,  std::vector<double> dev, std::vector<double> grad)
{
	m_Poly=pMesh;
    pMesh->NbFaceLabel=0;
    sort(grad.begin(), grad.end());
    //double epsilon = grad[grad.size()*1/30];
    vector<double>::iterator minimum;
    minimum = min_element (grad.begin(), grad.end());

    double threshold =  MedianofVector(grad);
    double thresholdDepth= MedianofVector(dev);
    int n_pixels = dev.size();

//    sort(grad.begin(), grad.end());
//    double thresholdGradientDepth = grad[grad.size()*1/7];

	 for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	    {
			pFacet->LabelVSA=0;
		}

	////Creation of the initial proxies by random seed triangle picking


// trouver les facets seeds
	m_Table_Proxy.clear();
    bool maxLOC= true;
    int i= 1;
     for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
		{

                Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
                Halfedge_around_facet_circulator end = pHalfedge;
                maxLOC = true;
				CGAL_For_all(pHalfedge,end)
				{
					Facet_handle pfacet2 = pHalfedge->opposite()->facet();
					if ((pfacet2 != NULL) && (pFacet->depthfacet < pfacet2->depthfacet) )
					{
					 maxLOC = false;
					}
				}

                if((maxLOC == true) && (pFacet != NULL) )///this triangle is chosen
                    {

                        //a proxy is created
                        Proxy NewProxy;
                        NewProxy.Normal=pFacet->normal();
                        pFacet->LabelVSA=i;
                        NewProxy.Seed=pFacet;
                        NewProxy.Label= i;
                        NewProxy.facets.push_back(pFacet);
                        ///the proxy is added
                        m_Table_Proxy.push_back(NewProxy);
                        pMesh->NbFaceLabel=  pMesh->NbFaceLabel +1;
                        i++;

                    }


	  }


        m_NbProxy=pMesh->NbFaceLabel;
        cout << "pMesh->NbFaceLabel =" << pMesh->NbFaceLabel << "et m_Table_Proxy.size()" << m_Table_Proxy.size() << endl;

        // Flooding
        std::list<Facet_iterator> seedFacets;
        int label;
        double indiceDepthREg, depthVoisin, indiceDepthREg2, gradientdepthVoisin ;
        Facet_iterator facet;
         bool first;
         double STD, STD2;
        std::vector <double> indiceDepthREgList;
        std::vector<Facet_handle> faces;
        std::vector<Proxy> ::iterator it;
        for(int i=0;i<m_Table_Proxy.size();i++)//For each proxy we select triangles that could grow
		{
			m_Table_Proxy[i].TabAdj.clear();
			Facet_iterator f=m_Table_Proxy[i].Seed;
			//seedFacets.push_front(f);
			seedFacets= addFacetToOrderedList(f, seedFacets);
			label= f->LabelVSA;
			//cout << "label" << label << endl;
            first = true;

            cout << "i =" << i << endl;
			while(!seedFacets.empty())
			{
                facet= seedFacets.front();
                seedFacets.pop_front();
                indiceDepthREgList.push_back(facet->depthfacet);
                Halfedge_around_facet_circulator pHalfedge = facet->facet_begin();
                Halfedge_around_facet_circulator end = pHalfedge;
                cout << "un" <<endl;
				CGAL_For_all(pHalfedge,end)
				{
					Facet_handle facet2 = pHalfedge->opposite()->facet();
					//depthVoisin= facet2->depthfacet;
					depthVoisin= facet2->depthfacet;
					gradientdepthVoisin= facet2->gradientdepth;
					//&& (first || (abs(depthVoisin - indiceDepthREg )< threshold3))
                    cout << "gradientdepthVoisin" << gradientdepthVoisin << endl;

                    //STD= standard_deviation_of(depth_facets_region(m_Table_Proxy[i].facets));



                    //|| (abs(depthVoisin - indiceDepthREg )> threshold3 )
					if ( (first) ||(gradientdepthVoisin< threshold) && (facet2->LabelVSA ==0))
					{
					cout << "kkk" << endl;
                        facet2->LabelVSA = label;
                       // seedFacets.push_front(facet2);
                        seedFacets= addFacetToOrderedList(facet2, seedFacets);
                        m_Table_Proxy[i].facets.push_back(facet2);

					}
					else if ((gradientdepthVoisin< threshold) && (facet2->LabelVSA != 0) && (facet2->LabelVSA != label))  // fusionner des régions
					{

                        	cout << "facet2->LabelVSA" <<facet2->LabelVSA<< "et label ="<< label <<endl;
////                  indice de la région du label facet2->LabelVSA
                      int ii = indiceRegion(facet2->LabelVSA);
                              cout << "ii =" << ii << "and i =" << i << endl;
                     //indiceDepthREg= MedianofVector(depth_facets_region(m_Table_Proxy[i].facets));
                    // indiceDepthREg2= MedianofVector(depth_facets_region(m_Table_Proxy[ii].facets));
                     //STD2= standard_deviation_of(depth_facets_region(m_Table_Proxy[ii].facets));
                    // Pour fusionner les deux régions, il faut ajouter les facettes de la région du label facet2->LabelVSA dans la région en cours et changer leurs label
                   if (ii != -1)
                   {
                  // cout << "le rapport =" << pow(mean_depth_facets_region(m_Table_Proxy[i].facets),2)/pow(mean_depth_facets_region(m_Table_Proxy[ii].facets),2) << endl;
                  double rapport = mean_depth_facets_region(m_Table_Proxy[i].facets)/mean_depth_facets_region(m_Table_Proxy[ii].facets);
                   //&& (abs(indiceDepthREg - indiceDepthREg2) > STD  || abs(indiceDepthREg - indiceDepthREg2) > STD2 )


                   //if ( (rapport > 0.75) && (rapport <1.25) (dR<dev) |||| (abs(indiceDepthREg - indiceDepthREg2) > STD  || abs(indiceDepthREg - indiceDepthREg2) > STD2)
                   if ((rapport > 0.50) && (rapport <1.75) )
                      {

                          // cout << "fusionnnnnnnnnnnnnnnnnnnnnnnnnn" << endl;

                           AddFacetsOfLabel1ToLabel2(ii, i);

                            it= m_Table_Proxy.begin() + ii;
                            // Supprimer cette région de m_NbProxy
                            it= m_Table_Proxy.erase(it);

                           //Supprimer tous les seeds qui ont label egal a facet2->LabelVSA
                            seedFacets = DeleteSeedRegionOfLabel(seedFacets, facet2->LabelVSA);

                            pMesh->NbFaceLabel=  pMesh->NbFaceLabel -1;
                            m_NbProxy = m_NbProxy-1;

                            if (ii < i)
                                i = i - 1;
                    }
                    }


                 }


				}
				first = false;

			}



		}
        cout << "first pre-processing" << endl;
		ReorganizeListRegions (pMesh);

        //  Filtering small regions
            for(int i=0;i<m_Table_Proxy.size();i++)//For each proxy we select triangles that could grow
		       {
                    faces= m_Table_Proxy[i].facets;

                    if (faces.size()< 3)
                    {

                        for (int j=0; j< faces.size(); j++)
                        {
                            faces.at(j)->LabelVSA= 0;
                        }
                                    it= m_Table_Proxy.begin() + i;
                                    // Supprimer cette région de m_NbProxy
                                    it= m_Table_Proxy.erase(it);
                                    pMesh->NbFaceLabel=  pMesh->NbFaceLabel-1;
                                    i=i-1;
                    }
                }

                // Fill holes
                cout << "m_Table_Proxy.size()" << m_Table_Proxy.size() << endl;


 //ConstructFaceColorMap (pMesh);
return pMesh;
}

PolyhedronPtr VSA_Component::Morpho_Fill_holes(PolyhedronPtr pMesh)
{
    double R, R2;
    std::vector<Facet_handle> facets;
    Point barycenterR, barycenterF ;
    double dist;

    for(int i=0;i<m_Table_Proxy.size();i++)//For each proxy we select triangles that could grow
     {
      cout << "region to fill " << i<< endl;

       R= std::sqrt(m_Table_Proxy[i].Area);
       R2= 2*R;
       pMesh= Morpho_Fill_holes_Using_Sphere(pMesh, i, R, R2);
       cout << "done " << endl;
       facets= m_Table_Proxy[i].facets;
       compute_region_center(pMesh,facets, barycenterR);
       for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
       {
         pMesh->compute_facet_center(pFacet,barycenterF);
         Vector deviation= barycenterR - barycenterF;
         dist= sqrt( pow(deviation.x(),2) + pow(deviation.y(),2) + pow(deviation.z(),2));
         if  (pFacet->morphofilled ==0 && (dist< R2))
         pFacet->LabelVSA=i;
       }
     }
     ReorganizeListRegions(pMesh);
     return pMesh;

}

PolyhedronPtr VSA_Component::Mathematical_morphology_post_processing(PolyhedronPtr pMesh)
{


 double area, RadiusVal;
 for(int i=0;i<m_Table_Proxy.size();i++)//For each proxy we select triangles that could grow
	{
//        if (i==10)
//        {
             area= m_Table_Proxy[i].Area;
             RadiusVal= std::sqrt(area);
             closing(pMesh,i, 0.1*RadiusVal);
             for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
            {
                if  (pFacet->morpho_closed !=0)
                     pFacet->LabelVSA=i;
            }



        //}


     }
//      for(int i=0;i<m_Table_Proxy.size();i++)//For each proxy we select triangles that could grow
//	{
////        if (i==10)
////        {
//
//                Fill_holes(pMesh,i);
//
//               for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
//            {
//                if  (pFacet->filled ==1)
//                     pFacet->LabelVSA=i;
//            }
//
//
//        //}
//
//
//     }

    ConstructFaceColorMap (pMesh);
return pMesh;
}

void VSA_Component::Burst_Wind_Segmentation_ON_Facets(PolyhedronPtr pMesh,  std::vector<double> dev)
{

	m_Poly=pMesh;
    pMesh->NbFaceLabel=0;
    sort(dev.begin(), dev.end());
    double epsilon = dev[dev.size()*1/5];
    double threshold =  mean_of(dev);

//    sort(grad.begin(), grad.end());
//    double thresholdGradientDepth = grad[grad.size()*1/7];

	 for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	    {
			pFacet->LabelVSA=0;
		}

	////Creation of the initial proxies by random seed triangle picking


// trouver les facets seeds
	m_Table_Proxy.clear();
    bool maxLOC= true;
    int i= 1;
     for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
		{

                Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
                Halfedge_around_facet_circulator end = pHalfedge;
                maxLOC = true;
				CGAL_For_all(pHalfedge,end)
				{
					Facet_handle pfacet2 = pHalfedge->opposite()->facet();
					if ((pfacet2 != NULL) && (pFacet->depthfacet < pfacet2->depthfacet) )
					{
					 maxLOC = false;
					}
				}

                if((maxLOC == true) && (pFacet != NULL) )///this triangle is chosen
                    {

                        //a proxy is created
                        Proxy NewProxy;
                        NewProxy.Normal=pFacet->normal();
                        pFacet->LabelVSA=i;
                        NewProxy.Seed=pFacet;
                        NewProxy.Label= i;
                        NewProxy.facets.push_back(pFacet);
                        ///the proxy is added
                        m_Table_Proxy.push_back(NewProxy);
                        pMesh->NbFaceLabel=  pMesh->NbFaceLabel +1;
                        i++;

                    }


	  }


        m_NbProxy=pMesh->NbFaceLabel;
        cout << "pMesh->NbFaceLabel =" << pMesh->NbFaceLabel << "et m_Table_Proxy.size()" << m_Table_Proxy.size() << endl;

        // Flooding
        std::list<Facet_iterator> seedFacets;
        int label;
        double indiceDepthREg, depthVoisin, indiceDepthREg2;
        Facet_iterator facet;
        double threshold3=0.0003;
         bool first;
         double STD;
        std::vector <double> indiceDepthREgList;
        for(int i=0;i<m_Table_Proxy.size();i++)//For each proxy we select triangles that could grow
		{
			m_Table_Proxy[i].TabAdj.clear();
			Facet_iterator f=m_Table_Proxy[i].Seed;
			seedFacets.push_front(f);
			label= f->LabelVSA;
			//cout << "label" << label << endl;
            first = true;


			while(!seedFacets.empty())
			{
                facet= seedFacets.front();
                seedFacets.pop_front();
                indiceDepthREgList.push_back(facet->depthfacet);
                Halfedge_around_facet_circulator pHalfedge = facet->facet_begin();
                Halfedge_around_facet_circulator end = pHalfedge;

				CGAL_For_all(pHalfedge,end)
				{
					Facet_handle facet2 = pHalfedge->opposite()->facet();
					depthVoisin= facet2->depthfacet;
					//&& (first || (abs(depthVoisin - indiceDepthREg )< threshold3))


                    STD= standard_deviation_of(depth_facets_region(m_Table_Proxy[i].facets));



                    //|| (abs(depthVoisin - indiceDepthREg )> threshold3 )
					if ( (first) ||(depthVoisin> threshold) && (facet2->LabelVSA ==0))
					{
                        facet2->LabelVSA = label;
                        seedFacets.push_front(facet2);
                        m_Table_Proxy[i].facets.push_back(facet2);

					}
					else if ((depthVoisin > threshold) && (facet2->LabelVSA != 0) && (facet2->LabelVSA != label))  // fusionner des régions
					{

//                   indice de la région du label facet2->LabelVSA
                     int ii = indiceRegion(facet2->LabelVSA);
                     indiceDepthREg= MedianofVector(depth_facets_region(m_Table_Proxy[i].facets));
                     indiceDepthREg2= MedianofVector(depth_facets_region(m_Table_Proxy[ii].facets));

                     //cout << "region en cours " << i << "sera fusionne avec" << ii << endl;
                    // Pour fusionner les deux régions, il faut ajouter les facettes de la région du label facet2->LabelVSA dans la région en cours et changer leurs label
                  //  cout << "before"<<m_Table_Proxy[i].facets.size() << endl;
                  cout << "le rapport =" << pow(mean_depth_facets_region(m_Table_Proxy[i].facets),2)/pow(mean_depth_facets_region(m_Table_Proxy[ii].facets),2) << endl;
                  if ( pow(mean_depth_facets_region(m_Table_Proxy[i].facets),2)/pow(mean_depth_facets_region(m_Table_Proxy[ii].facets),2)< 2)
                  {

                       AddFacetsOfLabel1ToLabel2(ii, i);

                         std::vector<Proxy> ::iterator it= m_Table_Proxy.begin() + ii;
                        // Supprimer cette région de m_NbProxy
                        it= m_Table_Proxy.erase(it);

                       //Supprimer tous les seeds qui ont label egal a facet2->LabelVSA
                        seedFacets = DeleteSeedRegionOfLabel(seedFacets, facet2->LabelVSA);

                        pMesh->NbFaceLabel=  pMesh->NbFaceLabel -1;
                        m_NbProxy = m_NbProxy-1;

                        if (ii < i)
                            i = i - 1;
                  }
                 }


				}
				first = false;

			}


		}

       ConstructFaceColorMap (pMesh);
}

double VSA_Component::variance_of(const std::vector<double>& values)
{
  double mean = mean_of(values);
  double sum = 0.0;
  for(unsigned long int i = 0; i < values.size(); ++i)
    sum += (values[i] - mean) * (values[i] - mean);
  return sum / static_cast<double>(values.size() - 1);
}



double VSA_Component::somme_carre(const std::vector<double>& values)
{
     double sum = 0.0;
  for(unsigned long int i = 0; i < values.size(); ++i)
  {
  sum += pow (values[i],2);
  }
  return sum;
}
// Fonction calculant l'écart-type des nombres contenus dans un vecteur
double VSA_Component::standard_deviation_of(const std::vector<double>& values)
{
  return sqrt(variance_of(values));
}


double VSA_Component::mean_depth_facets_region(std::vector<Facet_handle> facets)
{
     double nbr=0, sum=0, m;
     for (int i=0; i<facets.size(); i++)
     {
     nbr= nbr+1;
     sum=sum+ facets.at(i)->depthfacet;
     }
     return sum/nbr;
}

std::vector<double> VSA_Component::depth_facets_region(std::vector<Facet_handle> facets)
{
     std::vector<double> dep;
     for (int i=0; i<facets.size(); i++)
     {
        dep.push_back(facets.at(i)->depthfacet) ;
     }
     return dep;
}

 std::list<Facet_iterator> VSA_Component::DeleteSeedRegionOfLabel( std::list<Facet_iterator> seedFacets, int label )
{

    std::list<Facet_iterator> ::iterator it;
    Facet_iterator v;
    it=seedFacets.begin();
    while (it!=seedFacets.end())
    {
        v =*it;
         if (v->LabelVSA == label)
         it = seedFacets.erase(it);
         else ++it;
    }
    return seedFacets;

}

int VSA_Component::indiceRegion(int label)
{
  bool found= false;
   int i=0;

   while (!found && i<m_Table_Proxy.size())
   {
    if (m_Table_Proxy[i].Label == label)
    {
        found =true;

    }
     if (!found) i++;
   }

   if (m_Table_Proxy[i].Label == label)
    return i;
   else
    return -1;
}

void VSA_Component::AddFacetsOfLabel1ToLabel2(int label1, int label2)
{
   int i=0;

   std::vector<Facet_handle> f= m_Table_Proxy[label1].facets;

    assert(!f.empty());
   for (int j=0; j< f.size(); j++)
    {
     //cout << "here m_Table_Proxy[label2].Label" <<m_Table_Proxy[label2].Label << endl;
    // cout << "f.at(j)->LabelVSA" << f.at(j)->LabelVSA << endl;
         f.at(j)->LabelVSA = m_Table_Proxy[label2].Label;
        //  cout << "here" << endl;
    }

        if (f.size() > 1)
            m_Table_Proxy[label2].facets.insert(m_Table_Proxy[label2].facets.end(), f.begin(), f.end());
        else if (f.size() == 1)
        m_Table_Proxy[label2].facets.push_back(f.at(0));

}

double VSA_Component::mean_of(const std::vector<double>& values)
{
  double sum = 0.0;
  for(unsigned long int i = 0; i < values.size(); ++i)
    sum += values[i];
  return sum / static_cast<double>(values.size());
}

PolyhedronPtr VSA_Component::AdatedVersionOfTaubinSmoothing (PolyhedronPtr pMesh, double deformFactor, int iteraNum, bool preserveBoundaries)
{

Vertex_iterator	pVertex;
	int numVertex = pMesh->size_of_vertices();
	Vector * newPositions = new Vector[numVertex];

	for (int i=0; i<iteraNum; i++)
	{
		int n = 0;
		for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
		{

			Vector currentVector = pVertex->point() - CGAL::ORIGIN;
			Vector V_normal = pVertex->normal();

			// do not smooth the boundary vertices if demanded by user
			bool is_border_vertex = false;
			bool stopFlag = false;
			Halfedge_around_vertex_circulator hav = (*pVertex).vertex_begin();
			do
			{
				if (hav->is_border()==true)
				{
					is_border_vertex = true;
					stopFlag = true;
				}
				hav++;
			} while ((hav!=(*pVertex).vertex_begin())&&(stopFlag==false));

			if ((preserveBoundaries==true)&&(is_border_vertex==true))
			{
				newPositions[n] = currentVector;
				n++;
				continue;
			}

			std::size_t degree = (*pVertex).vertex_degree();
			double alpha = 1.0/degree;
			Vector vectemp = Point3d(0,0,0) - CGAL::ORIGIN;
			Vector vectNew;
			Halfedge_around_vertex_circulator h = (*pVertex).vertex_begin();

			do
			{
				vectemp = vectemp+(h->opposite()->vertex()->point()-CGAL::ORIGIN-currentVector)*alpha;
				++h;
			} while (h != (*pVertex).vertex_begin());


			double angle= Get_Angle_Weight (vectemp,V_normal);
            vectNew = currentVector + deformFactor*vectemp;
			if (angle < 0)
			{


                newPositions[n] = vectNew;
			}
			else
			{
                newPositions[n] = currentVector;
			}


			n++;
		}

		n = 0;
		for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
		{
			pVertex->point() = Point3d(0,0,0) + newPositions[n];
			n++;
		}
	}

	delete [] newPositions;
	newPositions = 0;

	pMesh->compute_normals();
	return pMesh;

}

PolyhedronPtr VSA_Component::LaplacianSmoothing (PolyhedronPtr pMesh, double deformFactor, int iteraNum, bool preserveBoundaries)
{

Vertex_iterator	pVertex;
	int numVertex = pMesh->size_of_vertices();
	Vector * newPositions = new Vector[numVertex];

	for (int i=0; i<iteraNum; i++)
	{
		int n = 0;
		for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
		{

			Vector currentVector = pVertex->point() - CGAL::ORIGIN;
			Vector V_normal = pVertex->normal();

			// do not smooth the boundary vertices if demanded by user
			bool is_border_vertex = false;
			bool stopFlag = false;
			Halfedge_around_vertex_circulator hav = (*pVertex).vertex_begin();
			do
			{
				if (hav->is_border()==true)
				{
					is_border_vertex = true;
					stopFlag = true;
				}
				hav++;
			} while ((hav!=(*pVertex).vertex_begin())&&(stopFlag==false));

			if ((preserveBoundaries==true)&&(is_border_vertex==true))
			{
				newPositions[n] = currentVector;
				n++;
				continue;
			}

			std::size_t degree = (*pVertex).vertex_degree();
			double alpha = 1.0/degree;
			Vector vectemp = Point3d(0,0,0) - CGAL::ORIGIN;
			Vector vectNew;
			Halfedge_around_vertex_circulator h = (*pVertex).vertex_begin();

			do
			{
				vectemp = vectemp+(h->opposite()->vertex()->point()-CGAL::ORIGIN-currentVector)*alpha;
				++h;
			} while (h != (*pVertex).vertex_begin());


			//double angle= Get_Angle_Weight (vectemp,V_normal);
            newPositions[n] = currentVector + deformFactor*vectemp;

			n++;
		}

		n = 0;
		for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
		{
			pVertex->point() = Point3d(0,0,0) + newPositions[n];
			n++;
		}
	}

	delete [] newPositions;
	newPositions = 0;

	pMesh->compute_normals();
	return pMesh;

}

PolyhedronPtr VSA_Component::Multi_scaleLaplacianSmoothing (PolyhedronPtr pMesh, double deformFactor, int iteraNum, bool preserveBoundaries, double radius)
{

    Vertex_iterator	pVertex, ppVertex;
	int numVertex = pMesh->size_of_vertices();
	Vector * newPositions = new Vector[numVertex];



    std::list<Vertex_handle> Vertexs_voisines;
    Vertex_handle vertex;

	for (int i=0; i<iteraNum; i++)
	{
		int n = 0;
		for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
		{
             // initialize to the vertices to visited ==0
             for (ppVertex = pMesh->vertices_begin(); ppVertex != pMesh->vertices_end(); ppVertex++)
              {
                ppVertex->VLS =0;
              }



			Vector currentVector = pVertex->point() - CGAL::ORIGIN;
			Vector V_normal = pVertex->normal();

			// do not smooth the boundary vertices if demanded by user
			bool is_border_vertex = false;
			bool stopFlag = false;
			Halfedge_around_vertex_circulator hav = (*pVertex).vertex_begin();
			do
			{
				if (hav->is_border()==true)
				{
					is_border_vertex = true;
					stopFlag = true;
				}
				hav++;
			} while ((hav!=(*pVertex).vertex_begin())&&(stopFlag==false));

			if ((preserveBoundaries==true)&&(is_border_vertex==true))
			{
				newPositions[n] = currentVector;
				n++;
				continue;
			}

			std::size_t degree = (*pVertex).vertex_degree();

			double cumul = 0;
			Vector vectemp = Point3d(0,0,0) - CGAL::ORIGIN;
			Vector vectNew;

			// considere the one ring
            Halfedge_around_vertex_circulator h = (*pVertex).vertex_begin();

            pVertex->VLS =1;
			do
			{
				double area = Area1Ring(h->opposite()->vertex());// aire de la 1-ring de h->opposite()->vertex()
				vectemp = vectemp+(h->opposite()->vertex()->point()-CGAL::ORIGIN-currentVector) * area;
				cumul += area;
				Vertexs_voisines.push_back(h->opposite()->vertex());
				h->opposite()->vertex()->VLS =1;
				++h;
			} while (h != (*pVertex).vertex_begin());

            // search for all the vertices neighbors of the current vertex inside a circle of radius r
            while(! Vertexs_voisines.empty() )
            {

             vertex= Vertexs_voisines.front();
             Vertexs_voisines.pop_front();
             Halfedge_around_vertex_circulator h = (*vertex).vertex_begin();
                 do
                {
                   Vector NVector =  h->opposite()->vertex()->point() - CGAL::ORIGIN;
                   Vector deviation= currentVector - NVector;
                   double dist= sqrt( deviation.x() * deviation.x() +
                                    deviation.y() * deviation.y() +
                                    deviation.z() * deviation.z());
                   //cout << "dist" << dist << endl;
                   if ((dist < radius) && (h->opposite()->vertex()->VLS == 0))
                   {
                    double area = Area1Ring(h->opposite()->vertex());// aire de la 1-ring de h->opposite()->vertex()
                    vectemp = vectemp+(h->opposite()->vertex()->point()-CGAL::ORIGIN-currentVector) * area;
                    cumul += area;
                    Vertexs_voisines.push_back(h->opposite()->vertex());
                   }
                    h->opposite()->vertex()->VLS =1;
                   ++h;
                } while (h != (*vertex).vertex_begin());
            }

            if (cumul != 0)
                vectemp = vectemp * 1. / cumul;

			double angle= Get_Angle_Weight (vectemp,V_normal);
			if (angle < 0)
			{
			 newPositions[n] = currentVector + deformFactor*vectemp;
			}
			else
			{
			  newPositions[n] = currentVector;
			}


			n++;
		}

		n = 0;
		for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
		{
			pVertex->point() = Point3d(0,0,0) + newPositions[n];
			n++;
		}
	}

	delete [] newPositions;
	newPositions = 0;

	pMesh->compute_normals();
	return pMesh;

}

double VSA_Component::Area1Ring(Vertex_handle pVertex)
{


           Halfedge_around_vertex_circulator h = (*pVertex).vertex_begin();
           double areas=0.0;
           Facet_handle pFacet1;
			do
			{
			    pFacet1 = h->opposite()->facet();
                double area = AreaFacetTriangle(pFacet1);
				areas += area;
				++h;
			} while (h != (*pVertex).vertex_begin());

         return areas;
}

double VSA_Component::Get_Angle_Weight(Vector V1 , Vector V2 )
{


	double Denom = std::sqrt(V1*V1) * std::sqrt(V2*V2);

	double Multi = V1*V2;

	double Cosine = 0;
	if (Denom != 0.0)
		Cosine = Multi / Denom;
	else
		return 0.0;

	if (Cosine > 1.0)
		Cosine = 1.0;
	if (Cosine < -1.0)
		Cosine = -1.0;

	return Cosine;
}

 PolyhedronPtr VSA_Component::AdjustDepth(PolyhedronPtr pMesh)
 {

int n= 0;
double mindepth;
Vertex_iterator	pVertex;
  //Find the smallest depth(non null)
  for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
		{
            if (n==0)
            {
                mindepth= pVertex->depthvertex;
            }
            else
            {
              if (mindepth > pVertex->depthvertex)
              {
               mindepth= pVertex->depthvertex;
              }
            }
            n++;
		}

		int i=1;
		for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
		{
            pVertex->depthvertex= pVertex->depthvertex + mindepth * (i/n);
            i++;
		}


    return pMesh;
 }

std::vector<double> VSA_Component::computeGradientMagnitudeFacet(PolyhedronPtr* pMesh)
{


   std::vector<double> grad;
    for (Facet_iterator pFacet = (*pMesh)->facets_begin(); pFacet != (*pMesh)->facets_end(); pFacet++)
	{
        Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
		Point3d P = pHalfedge->vertex()->point();
		Point3d Q = pHalfedge->next()->vertex()->point();
		Point3d R = pHalfedge->next()->next()->vertex()->point();

		Vector deviation= P - Q;
		double ax= deviation.x();
		double ay= deviation.y();
		double az= deviation.z();

		// change of the depth between P and Q
		double da= (pHalfedge->vertex()->depthvertex) - (pHalfedge->next()->vertex()->depthvertex);


		Vector deviation2= P - R;
		double bx= deviation2.x();
		double by= deviation2.y();
		double bz= deviation2.z();

		// change of the depth between P and R
        double db= (pHalfedge->vertex()->depthvertex) - (pHalfedge->next()->next()->vertex()->depthvertex);


		// Normal de la facette
		Vector normal = pFacet->normal();
		double nx= normal.x();
		double ny= normal.y();
		double nz= normal.z();

        Matrix A(3,3);
        Vector_OpenNL X(3);
        Vector_OpenNL B(3);

         A.set_coef(0,0,ax);
         A.set_coef(0,1,ay);
         A.set_coef(0,2,az);

         A.set_coef(1,0,bx);
         A.set_coef(1,1,by);
         A.set_coef(1,2,bz);

         A.set_coef(2,0,nx);
         A.set_coef(2,1,ny);
         A.set_coef(2,2,nz);
         B[0]=da;
         B[1]=db;
         B[2]=0.0;
         X[0] = X[1] = X[2] = 0.0;

        double d;

       SparseLA m_linearAlgebra;

       m_linearAlgebra.linear_solver(A,B,X,d);

       pFacet->gradientdepth = sqrt( pow(X[0],2) + pow(X[1],2) + pow(X[2],2))  ;
       if (pFacet->gradientdepth <0 || isnan(pFacet->gradientdepth))
       {
       cout << "waaak gradient" << endl;
       }
      grad.push_back(pFacet->gradientdepth);

	}
	return grad;

}



VSA_Component::VSA_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p)
{
	int i=0;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.515600;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.531300;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.546900;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.562500;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.578100;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.593800;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.609400;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.625000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.640600;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.656300;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.671900;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.687500;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.703100;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.718800;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.734400;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.750000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.765600;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.781300;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.796900;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.812500;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.828100;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.843800;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.859400;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.875000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.890600;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.906300;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.921900;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.937500;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.953100;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.968800;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.984400;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	1.000000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.015600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.031300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.046900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.062500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.078100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.093800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.109400;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.125000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.140600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.156300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.171900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.187500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.203100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.218800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.234400;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.250000;	LUT_Seg[i++]=	1.000000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.265600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.281300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.296900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.312500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.328100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.343800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.359400;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.375000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.390600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.406300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.421900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.437500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.453100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.468800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.484400;	LUT_Seg[i++]=	1.000000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.500000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.515600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.531300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.546900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.562500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.578100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.593800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.609400;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.625000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.640600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.656300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.671900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.687500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.703100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.718800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.734400;	LUT_Seg[i++]=	1.000000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.750000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.765600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.781300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.796900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.812500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.828100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.843800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.859400;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.875000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.890600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.906300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.921900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.937500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.953100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.968800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.984400;	LUT_Seg[i++]=	1.000000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.015600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.031300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.984400;		LUT_Seg[i++]=	0.046900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.968800;		LUT_Seg[i++]=	0.062500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.953100;		LUT_Seg[i++]=	0.078100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.937500;		LUT_Seg[i++]=	0.093800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.921900;		LUT_Seg[i++]=	0.109400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.906300;		LUT_Seg[i++]=	0.125000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.890600;		LUT_Seg[i++]=	0.140600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.875000;		LUT_Seg[i++]=	0.156300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.859400;		LUT_Seg[i++]=	0.171900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.843800;		LUT_Seg[i++]=	0.187500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.828100;		LUT_Seg[i++]=	0.203100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.812500;		LUT_Seg[i++]=	0.218800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.796900;		LUT_Seg[i++]=	0.234400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.781300;
	LUT_Seg[i++]=	0.250000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.765600;		LUT_Seg[i++]=	0.265600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.750000;		LUT_Seg[i++]=	0.281300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.734400;		LUT_Seg[i++]=	0.296900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.718800;		LUT_Seg[i++]=	0.312500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.703100;		LUT_Seg[i++]=	0.328100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.687500;		LUT_Seg[i++]=	0.343800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.671900;		LUT_Seg[i++]=	0.359400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.656300;		LUT_Seg[i++]=	0.375000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.640600;		LUT_Seg[i++]=	0.390600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.625000;		LUT_Seg[i++]=	0.406300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.609400;		LUT_Seg[i++]=	0.421900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.593800;		LUT_Seg[i++]=	0.437500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.578100;		LUT_Seg[i++]=	0.453100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.562500;		LUT_Seg[i++]=	0.468800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.546900;		LUT_Seg[i++]=	0.484400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.531300;
	LUT_Seg[i++]=	0.500000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.515600;		LUT_Seg[i++]=	0.515600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.500000;		LUT_Seg[i++]=	0.531300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.484400;		LUT_Seg[i++]=	0.546900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.468800;		LUT_Seg[i++]=	0.562500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.453100;		LUT_Seg[i++]=	0.578100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.437500;		LUT_Seg[i++]=	0.593800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.421900;		LUT_Seg[i++]=	0.609400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.406300;		LUT_Seg[i++]=	0.625000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.390600;		LUT_Seg[i++]=	0.640600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.375000;		LUT_Seg[i++]=	0.656300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.359400;		LUT_Seg[i++]=	0.671900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.343800;		LUT_Seg[i++]=	0.687500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.328100;		LUT_Seg[i++]=	0.703100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.312500;		LUT_Seg[i++]=	0.718800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.296900;		LUT_Seg[i++]=	0.734400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.281300;
	LUT_Seg[i++]=	0.750000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.265600;		LUT_Seg[i++]=	0.765600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.250000;		LUT_Seg[i++]=	0.781300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.234400;		LUT_Seg[i++]=	0.796900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.218800;		LUT_Seg[i++]=	0.812500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.203100;		LUT_Seg[i++]=	0.828100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.187500;		LUT_Seg[i++]=	0.843800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.171900;		LUT_Seg[i++]=	0.859400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.156300;		LUT_Seg[i++]=	0.875000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.140600;		LUT_Seg[i++]=	0.890600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.125000;		LUT_Seg[i++]=	0.906300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.109400;		LUT_Seg[i++]=	0.921900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.093800;		LUT_Seg[i++]=	0.937500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.078100;		LUT_Seg[i++]=	0.953100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.062500;		LUT_Seg[i++]=	0.968800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.046900;		LUT_Seg[i++]=	0.984400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.031300;
	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.015600;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.984400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.968800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.953100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.937500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.921900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.906300;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.890600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.875000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.859400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.843800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.828100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.812500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.796900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.781300;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.765600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.750000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.734400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.718800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.703100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.687500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.671900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.656300;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.640600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.625000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.609400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.593800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.578100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.562500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.546900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.531300;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.515600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.500000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.484400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.468800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.453100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.437500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.421900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.406300;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.390600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.375000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.359400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.343800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.328100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.312500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.296900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.281300;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.265600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.250000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.234400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.218800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.203100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.187500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.171900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.156300;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.140600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.125000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.109400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.093800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.078100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.062500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.046900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.031300;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.015600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.984400;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.968800;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.953100;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.937500;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.921900;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.906300;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.890600;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.875000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.859400;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.843800;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.828100;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.812500;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.796900;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.781300;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	0.765600;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.750000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.734400;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.718800;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.703100;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.687500;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.671900;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.656300;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.640600;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.625000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.609400;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.593800;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.578100;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.562500;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.546900;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.531300;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	0.515600;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;


	// from IHM
	displayFaceLabel=false;

	// MEPP 2
	componentName = "VSA_Component";
	init = 1;
}

void VSA_Component::Flooding()
	{
		m_Poly->NbFaceLabel=m_NbProxy;
		typedef std::multiset<FacetToIntegrate,CompFacet> ListFacet_model;
		typedef std::multiset<FacetToIntegrate,CompFacet>::iterator ListFacet_iterator;

		for(Facet_iterator	pface	=	m_Poly->facets_begin();
				pface	!= m_Poly->facets_end();
				pface++)
		{
			pface->LabelVSA=-1;
		}


		ListFacet_model ListFacet;
		for(int i=0;i<m_NbProxy;i++)//For each proxy we select triangles that could grow
		{
			m_Table_Proxy[i].TabAdj.clear();
			Facet_iterator f=m_Table_Proxy[i].Seed;
			f->LabelVSA=i;


			//we extract the three triangles
			Facet_iterator ff1, ff2, ff3;
			FacetToIntegrate f1, f2, f3;

			Halfedge_around_facet_circulator pHalfedge = f->facet_begin();

			if(!pHalfedge->opposite()->is_border())
			{
				ff1 = pHalfedge->opposite()->facet();
				f1.Facet=ff1;
				f1.PossibleCluster=i;
				f1.DistanceLLoyd=DistorsionError(f1.Facet,m_Table_Proxy[i]);
				ListFacet.insert(f1);

			}
			if(!pHalfedge->next()->opposite()->is_border())
			{
				ff2 = pHalfedge->next()->opposite()->facet();
				f2.Facet=ff2;
				f2.PossibleCluster=i;
				f2.DistanceLLoyd=DistorsionError(f2.Facet,m_Table_Proxy[i]);
				ListFacet.insert(f2);
			}


			if(!pHalfedge->next()->next()->opposite()->is_border())
			{
				ff3 = pHalfedge->next()->next()->opposite()->facet();
				f3.Facet=ff3;
				f3.PossibleCluster=i;
				f3.DistanceLLoyd=DistorsionError(f3.Facet,m_Table_Proxy[i]);
				ListFacet.insert(f3);
			}
		}

		ListFacet_iterator it;
		for(it=ListFacet.begin();it!=ListFacet.end();)
		{
			if(it->Facet->LabelVSA==-1)
			{
				it->Facet->LabelVSA=it->PossibleCluster;
				//we add adjacent triangles to the queue


				//we extract the three triangles
				Facet_iterator ff1, ff2, ff3;

				Halfedge_around_facet_circulator pHalfedge = it->Facet->facet_begin();


				FacetToIntegrate f1, f2, f3;
				if(!pHalfedge->opposite()->is_border())
				{
					ff1 = pHalfedge->opposite()->facet();
					if(ff1->LabelVSA==-1 )
					{
						f1.Facet=ff1;
						f1.PossibleCluster=it->PossibleCluster;
						f1.DistanceLLoyd=DistorsionError(f1.Facet,m_Table_Proxy[it->PossibleCluster]);
						ListFacet.insert(f1);
					}
				}
				if(!pHalfedge->next()->opposite()->is_border())
				{
					ff2 = pHalfedge->next()->opposite()->facet();
					if(ff2->LabelVSA==-1)
					{
						f2.Facet=ff2;
						f2.PossibleCluster=it->PossibleCluster;
						f2.DistanceLLoyd=DistorsionError(f2.Facet,m_Table_Proxy[it->PossibleCluster]);
						ListFacet.insert(f2);
					}
				}
				if(!pHalfedge->next()->next()->opposite()->is_border())
				{
					ff3 = pHalfedge->next()->next()->opposite()->facet();
					if(ff3->LabelVSA==-1)
					{
						f3.Facet=ff3;
						f3.PossibleCluster=it->PossibleCluster;
						f3.DistanceLLoyd=DistorsionError(f3.Facet,m_Table_Proxy[it->PossibleCluster]);
						ListFacet.insert(f3);
					}
				}
			}
			ListFacet.erase(it);
			it=ListFacet.begin();


		}

		//integration of the adjacency information between proxies
	/*	for(Halfedge_iterator pHalfedge	=	m_Poly->halfedges_begin();
				pHalfedge	!= m_Poly->halfedges_end();
				pHalfedge++)
		{
			if(pHalfedge->is_border()||pHalfedge->opposite()->is_border())
				continue;

			int Label1=pHalfedge->facet()->LabelVSA;
			int Label2=pHalfedge->opposite()->facet()->LabelVSA;
			if(Label1!=Label2)
			{

				bool IsFound=false;
				for(unsigned int i=0;i<m_Table_Proxy[Label1].TabAdj.size();i++)
					if(m_Table_Proxy[Label1].TabAdj[i]==Label2)
						IsFound=true;
				if(IsFound==false)
				{
					m_Table_Proxy[Label1].TabAdj.push_back(Label2);
					m_Table_Proxy[Label2].TabAdj.push_back(Label1);

				}


			}


		}*/



	}

	void VSA_Component::ProxyFitting()
	{
		Vector * TabNormal=new Vector[m_NbProxy];
		double * TabArea=new double[m_NbProxy];

		double * DistanceMin=new double[m_NbProxy];
		double * DistanceMax=new double[m_NbProxy];

		for (int i=0;i<m_NbProxy;i++)
		{
			TabArea[i]=0;
			TabNormal[i]=CGAL::NULL_VECTOR;
			DistanceMin[i]=100000000;
			DistanceMax[i]=0;



		}
		for(Facet_iterator	pface	=	m_Poly->facets_begin();
				pface	!= m_Poly->facets_end();
				pface++)
		{
			double area=AreaFacetTriangleSeg(pface);
			TabArea[pface->LabelVSA]+=area;
			TabNormal[pface->LabelVSA]=TabNormal[pface->LabelVSA]+pface->normal()*area;

		}

		for (int i=0;i<m_NbProxy;i++)
		{

			m_Table_Proxy[i].Normal=TabNormal[i]/TabArea[i];
			m_Table_Proxy[i].Area=TabArea[i];
			m_Table_Proxy[i].TotalDistorsion=0;

		}



		// a new seed is assigned to each proxy
		for(Facet_iterator	pface	=	m_Poly->facets_begin();
				pface	!= m_Poly->facets_end();
				pface++)
		{
			double distance=DistorsionError(pface,m_Table_Proxy[pface->LabelVSA]);
			m_Table_Proxy[pface->LabelVSA].TotalDistorsion+=distance;
			if(distance<DistanceMin[pface->LabelVSA])
			{

				DistanceMin[pface->LabelVSA]=distance;
				m_Table_Proxy[pface->LabelVSA].Seed=pface;
			}

			//we pick the facet corresponding to the max distorsion
			if(distance>DistanceMax[pface->LabelVSA])
			{

				DistanceMax[pface->LabelVSA]=distance;
				m_Table_Proxy[pface->LabelVSA].MostDistordedFacet=pface;
			}


		}

		delete []DistanceMin;
		delete []TabNormal;
		delete []TabArea;
		delete []DistanceMax;



	}


	void VSA_Component::ProxyInsertion()
	{
		int NumProxMax=0;
		double DistorsionMax=0;

		EvalInsertion(NumProxMax,DistorsionMax);
		CreateNewProxy( NumProxMax);

	}

	void VSA_Component::EvalInsertion(int & NumProxMax,double & DistorsionMax)
	{
		for (int i=0;i<m_NbProxy;i++)
		{

			if(	m_Table_Proxy[i].TotalDistorsion>DistorsionMax)
			{
				NumProxMax=i;
				DistorsionMax=m_Table_Proxy[i].TotalDistorsion;

			}
		}

	}

	void VSA_Component::CreateNewProxy(int NumProxMax)
	{
		Proxy NewProxy;
		if(m_Table_Proxy[NumProxMax].MostDistordedFacet!=m_Table_Proxy[NumProxMax].Seed)
			NewProxy.Seed=m_Table_Proxy[NumProxMax].MostDistordedFacet;
		else
		{
			m_Table_Proxy[NumProxMax].MostDistordedFacet++;
			NewProxy.Seed=m_Table_Proxy[NumProxMax].MostDistordedFacet;
		}


		NewProxy.Normal=m_Table_Proxy[NumProxMax].MostDistordedFacet->normal();
		m_NbProxy++;
		m_Table_Proxy.push_back(NewProxy);
		NewProxy.Seed->tag(15);
	}


	void VSA_Component::Variational_SegmentationIncr(PolyhedronPtr pMesh, int NbRegions, int NbIter)
	{
		m_Poly=pMesh;
		Init(2);
		Flooding();

		for(int i=0;i<NbRegions-2;i++)
		{
			for(int j=0;j<NbIter;j++)
			{
				ProxyFitting();
				Flooding();
			}
			ProxyFitting();
			ProxyInsertion();
			Flooding();

		}
		for(int j=0;j<NbIter;j++)
		{
			ProxyFitting();
			Flooding();
		}

	}

	void VSA_Component::Variational_Segmentation(PolyhedronPtr pMesh, int NbRegions, int NbIter)
	{
		m_Poly=pMesh;
		Init(NbRegions);
		Flooding();


		for(int i=0;i<NbIter;i++)
		{
			ProxyFitting();
			Flooding();
		}

	}

	void VSA_Component::ConstructFaceDepthMap(PolyhedronPtr pMesh)
	{

		Facet_iterator pFacet	=	pMesh->facets_begin();
        double Min = pFacet->depthfacet;
		double Max = pFacet->depthfacet;
		pFacet++;
		for(;pFacet	!= pMesh->facets_end();pFacet++)
		{
		 if(Min > pFacet->depthfacet) Min = pFacet->depthfacet;
         if(Max < pFacet->depthfacet) Max = pFacet->depthfacet;
		}
        cout << "Min " << Min << endl;
        cout << "Max" << Max << endl;
		pFacet	=	pMesh->facets_begin();
		for(;pFacet	!= pMesh->facets_end();pFacet++)
		{

			double R=(double)(pFacet->depthfacet)/(Max - Min)*255.;
			int indiceLut=floor(R);

			pFacet->color(LUT_Seg[3*indiceLut],LUT_Seg[3*indiceLut+1],LUT_Seg[3*indiceLut+2]);

		}
	}

		void VSA_Component::ConstructGradFaceDepthMap(PolyhedronPtr pMesh, double threshold )
	{

		Facet_iterator pFacet	=	pMesh->facets_begin();
        double Min = pFacet->gradientdepth;
		double Max = pFacet->gradientdepth;
		pFacet++;
		for(;pFacet	!= pMesh->facets_end();pFacet++)
		{
		 if(Min > pFacet->gradientdepth) Min = pFacet->gradientdepth;
         if(Max < pFacet->gradientdepth) Max = pFacet->gradientdepth;
		}
        cout << "Min " << Min << endl;
        cout << "Max" << Max << endl;
		pFacet	=	pMesh->facets_begin();
		for(;pFacet	!= pMesh->facets_end();pFacet++)
		{

//			double R=(double)(pFacet->gradientdepth)/(Max - Min)*255.;
//			int indiceLut=(int)floor(R)% 256;
//
//			pFacet->color(LUT_Seg[3*indiceLut],LUT_Seg[3*indiceLut+1],LUT_Seg[3*indiceLut+2]);
            if (pFacet->gradientdepth > threshold)
            pFacet->color(1,0,0);
            else
            pFacet->color(0,1,0);
		}
	}

	void VSA_Component::ConstructFaceColorMap(PolyhedronPtr pMesh)
{

                //Vertex_iterator pVertex = NULL; // MT

		Facet_iterator pFacet	=	pMesh->facets_begin();
		for(;pFacet	!= pMesh->facets_end();pFacet++)
		{

			double R=(double)(pFacet->LabelVSA)/(double)pMesh->NbFaceLabel*255.;
			int indiceLut=(int)floor(R) % 256;

            if (indiceLut == 0)
                pFacet->color(0, 0, 0);
            else
                pFacet->color(LUT_Seg[3*indiceLut],LUT_Seg[3*indiceLut+1],LUT_Seg[3*indiceLut+2]);

		}
}

	void VSA_Component::ConstructClosedFaceColorMap(PolyhedronPtr pMesh)
{

		Facet_iterator pFacet	=	pMesh->facets_begin();
		for(;pFacet	!= pMesh->facets_end();pFacet++)
		{
		if (pFacet->morpho_closed ==1)
					pFacet->color(1,0,0);
		else
					pFacet->color(0,1,0);



		}
                //Vertex_iterator pVertex = NULL; // MT

//		Facet_iterator pFacet	=	pMesh->facets_begin();
//		for(;pFacet	!= pMesh->facets_end();pFacet++)
//		{
//
//			double R=(double)(pFacet->closed)/(double)2.0*255.;
//			int indiceLut=(int)floor(R) % 256;
//
//			pFacet->color(LUT_Seg[3*indiceLut],LUT_Seg[3*indiceLut+1],LUT_Seg[3*indiceLut+2]);
//
//		}
}

void VSA_Component::ConstructErodedFaceColorMap(PolyhedronPtr pMesh)
{

                //Vertex_iterator pVertex = NULL; // MT

		Facet_iterator pFacet	=	pMesh->facets_begin();
		for(;pFacet	!= pMesh->facets_end();pFacet++)
		{
		if (pFacet->eroded ==1)
					pFacet->color(1,0,0);
		else
					pFacet->color(0,1,0);



		}
}

	void VSA_Component::ConstructDilatedFaceColorMap(PolyhedronPtr pMesh)
{

            Facet_iterator pFacet	=	pMesh->facets_begin();
            for(;pFacet	!= pMesh->facets_end();pFacet++)
            {
            if (pFacet->morpho_dilated ==1)
                        pFacet->color(1,0,0);
            else
                        pFacet->color(0,1,0);



            }
                //Vertex_iterator pVertex = NULL; // MT

//		Facet_iterator pFacet	=	pMesh->facets_begin();
//		for(;pFacet	!= pMesh->facets_end();pFacet++)
//		{
//
//			double R=(double)(pFacet->dilated)/(double)2.0*255.;
//			int indiceLut=(int)floor(R) % 256;
//
//			pFacet->color(LUT_Seg[3*indiceLut],LUT_Seg[3*indiceLut+1],LUT_Seg[3*indiceLut+2]);
//
//		}
}
#endif
