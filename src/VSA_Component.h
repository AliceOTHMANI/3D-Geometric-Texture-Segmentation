/*!
	\file VSA_Component.h
	\brief Segmentation of a polyhedra into planar proxies
	\brief According to: Variational Shape Approximation
		David Cohen-Steiner, Pierre Alliez and Mathieu Desbrun, SIGGRAPH '2004.

	\author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
	\date	2008
 */
#ifndef VSA_COMPONENT_H
#define VSA_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_VSA

#include "../../../../mepp/mepp_component.h"
#include "VSA_Polyhedron.h"

class Viewer;

/**
 \class	VSA_Component

 \brief	Segmentation of a polyhedra


 */
class VSA_Component :
  public mepp_component
{
	public:

		/**
		 \struct	Proxy

		 \brief	Structure of a proxy

		 */
		 typedef CGAL::Simple_cartesian<double>					AABB_Kernel;
         typedef AABB_Kernel::Point_3                            Point;
		typedef struct
		  {

			Vector Normal;///< The normal of the proxy

			Facet_iterator Seed;///< The seed of the proxy
			double Area; ///< The area of the proxy
			std::vector<int> TabAdj; ///< Indices of the neighbooring proxies
			double TotalDistorsion; ///< Global approximation error over the proxy
            std::vector<Facet_handle> facets; /// the facet of the region
			Facet_iterator MostDistordedFacet; ///< The most distorded facet
			int Label; ///< Proxy label
		  }
		  Proxy;



		/*!
		* \brief Constructor
		*/
		VSA_Component(Viewer* v, PolyhedronPtr p);

		/*!
		* \brief Destructor
		*/
		~VSA_Component() {}

		/**
		 \fn	void VSA_Component::Variational_Segmentation(PolyhedronPtr pMesh,int NbRegions,
		 		int NbIter);

		 \brief	Variational segmentation, classical version: seeds are sampled and then relaxed

		 \param	pMesh	 	The mesh.
		 \param	NbRegions	The number of regions.
		 \param	NbIter   	The number relaxation iterations.
		 */
		void Variational_Segmentation(PolyhedronPtr pMesh,int NbRegions, int NbIter);

		/**
		 \fn	void VSA_Component::Variational_SegmentationIncr(PolyhedronPtr pMesh,int NbRegions,
		 		int NbIter);

		 \brief	Variational segmentation, incremental version: seeds are incrementally added and relaxed.

		 Longer but more accurate than the classical version.

		 \param	pMesh	 	The mesh.
		 \param	NbRegions	The number of proxies.
		 \param	NbIter   	The number relaxation iterations.
		 */
		void Variational_SegmentationIncr(PolyhedronPtr pMesh,int NbRegions, int NbIter);//algorithme icrémental

		/**
		 \fn	void VSA_Component::Init(int NbProxy);

		 \brief	Initialises the process.

		 \param	NbProxy	The number of proxies.
		 */

		void Init(int NbProxy);

		/**
		 \fn	void VSA_Component::Flooding();

		 \brief	Flooding process. See the Variational Shape Approximation paper for details.

		 */

		void Flooding();

		/**
		 \fn	void VSA_Component::ProxyFitting();

		 \brief	Proxy fitting. Triangles are integrated in their respective proxy.

		 See the Variational Shape Approximation paper for details.



		 */

		void ProxyFitting();

		/**
		 \fn	double VSA_Component::DistorsionError(Facet_iterator f,Proxy p);

		 \brief	Evaluate the difference of normal between the facet and the proxy


		 \param	f	The facet
		 \param	p	The proxy

		 \return	.
		 */
		double DistorsionError(Facet_iterator f,Proxy p);

		/**
		 \fn	void VSA_Component::ProxyInsertion();

		 \brief	Insertion of a new proxy

		 */

		void ProxyInsertion();

		/**
		 \fn	void VSA_Component::EvalInsertion(int & NumProxMax,double & DistorsionMax);

		 \brief	Find the proxy where the approximation error is maximum, i.e. where an insertion of a new proxy is necessary


		 \param [out]	NumProxMax   	Label of the proxy associated with the highest error
		 \param [out]	DistorsionMax	Value of the maximum error
		 */

		void EvalInsertion(int & NumProxMax,double & DistorsionMax);

		/**
		 \fn	void VSA_Component::CreateNewProxy(int NumProxMax);

		 \brief	Creates a new proxy starting from the highest error facet from proxy NumProxMax

		\param	NumProxMax   	Label of the proxy associated with the highest error
		 */

		void CreateNewProxy(int NumProxMax);

		/**
		 \fn	void VSA_Component::ConstructFaceColorMap(PolyhedronPtr pMesh);

		 \brief	This method map the segmentation into facet colors

		 \param	pMesh	The mesh.
		 */

		void ConstructFaceColorMap(PolyhedronPtr pMesh);//creation de la table couleur pour l'affichage des régions

		/** The Burst wind segmentation methods
		*/
		PolyhedronPtr AdatedVersionOfTaubinSmoothing (PolyhedronPtr pMesh, double deformFactor, int iteraNum, bool preserveBoundaries);
		PolyhedronPtr LaplacianSmoothing (PolyhedronPtr pMesh, double deformFactor, int iteraNum, bool preserveBoundaries);
		PolyhedronPtr Multi_scaleLaplacianSmoothing (PolyhedronPtr pMesh, double deformFactor, int iteraNum, bool preserveBoundaries, double radius);
		double Get_Angle_Weight(Vector V1 , Vector V2 );
        double Area1Ring(Vertex_handle pVertex);
		double getMaxDim(PolyhedronPtr polyhedron_ptr);
		std::vector<double> computeLaplacianCotg2(PolyhedronPtr pMesh);
		std::vector<double> computeLaplacianCotg3(PolyhedronPtr pMesh);
        double computeWcot(Halfedge_handle h);
        double computeHalfWcot(Halfedge_handle h);
        double acos_stable(double CosAng);
        void computeGaussianCurv1Ring(PolyhedronPtr pMesh);
        double getAngle(Vector he1,Vector he2);
        std::vector<double> DeviationPolyhedronToPolyhedron(PolyhedronPtr original, PolyhedronPtr lisse);
        PolyhedronPtr DepthFacets(PolyhedronPtr pMesh);
        void ConstructFaceDepthMap(PolyhedronPtr pMesh);
        void WriteTxt2(std::vector<double> dev, std::string file);
        std::vector<double> computeGradientMagnitudeFacet(PolyhedronPtr* pMesh);
        void Burst_Wind_Segmentation_ON_Facets(PolyhedronPtr pMesh, std::vector<double> dev);
        PolyhedronPtr Burst_Wind_Segmentation_ON_Facets_2(PolyhedronPtr pMesh,  std::vector<double> dev, std::vector<double> grad);
        double MedianofVector (std::vector<double> grades);
        double mean_of(const std::vector<double>& values);
        void AddFacetsOfLabel1ToLabel2(int, int);
        int indiceRegion(int);
        std::list<Facet_iterator> DeleteSeedRegionOfLabel( std::list<Facet_iterator> Vertexs, int label );
        double mean_depth_facets_region(std::vector<Facet_handle> facets);
        double standard_deviation_of(const std::vector<double>& values);
        double variance_of(const std::vector<double>& values);
        std::vector<double> depth_facets_region(std::vector<Facet_handle> facets);
        double getOtsuThreshold(std::vector<double> data);
        std::list<Facet_iterator> addFacetToOrderedList(Facet_iterator pVertex,  std::list<Facet_iterator> Vertexs );
        double somme_carre(const std::vector<double>& values);
        PolyhedronPtr dilatation(PolyhedronPtr pMesh,int indR,double R );
        void Erosion(PolyhedronPtr pMesh, int indR, double R );
        PolyhedronPtr Erosion_after_dilation(PolyhedronPtr pMesh, int indR, double R );
        PolyhedronPtr closing(PolyhedronPtr pMesh, int indR, double R );
        PolyhedronPtr Morpho_Fill_holes(PolyhedronPtr pMesh, int indR);
        PolyhedronPtr Mathematical_morphology_post_processing(PolyhedronPtr pMesh);
        void ConstructDilatedFaceColorMap(PolyhedronPtr pMesh);
        void ReorganizeListRegions(PolyhedronPtr pMesh);
        double AreaFacetTriangle(Facet_handle &f);
        bool distFtoList(PolyhedronPtr pMesh, Facet_handle face, std::vector<Facet_handle> facets, double R);
        void ConstructErodedFaceColorMap(PolyhedronPtr pMesh);
        void ConstructClosedFaceColorMap(PolyhedronPtr pMesh);
        void ConstructGradFaceDepthMap(PolyhedronPtr pMesh, double threshold);
        PolyhedronPtr TranslateToZeroTheDepth(PolyhedronPtr pMesh, std::vector<double> dev);
        PolyhedronPtr AdjustDepth(PolyhedronPtr pMesh);
        void compute_region_center(PolyhedronPtr pMesh, std::vector<Facet_handle> facets, Point& center);
        PolyhedronPtr Morpho_Fill_holes_Using_Sphere(PolyhedronPtr pMesh, int indR, double R, double R2);
        PolyhedronPtr Morpho_Fill_holes(PolyhedronPtr pMesh);
        void Get_Clicked_Vertices(PolyhedronPtr pMesh, double x, double y, int tolerance);
        void draw_Sphere(PolyhedronPtr pMesh, float r, int sphereDisk);

	public:


		int m_NbProxy; ///< The number of proxies

		PolyhedronPtr m_Poly;///< The mesh

		std::vector<Proxy> m_Table_Proxy; ///< The set of proxies

		double LUT_Seg[3*256]; ///< Look-up table for color rendering

	// from IHM
	private:

		bool displayFaceLabel; ///< Boolean indicating activation or deactivation of the segmentation rendering
};

#endif

#endif // VSA_COMPONENT_H
