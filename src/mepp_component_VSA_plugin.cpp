#include <mepp_config.h>
#ifdef BUILD_component_VSA

#include "mepp_component_VSA_plugin.hxx"

#include "dialSettings.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>

#include "VSA_Component.h"
typedef boost::shared_ptr<VSA_Component> VSA_ComponentPtr;

#include "../../../Analysis/Curvature/src/Curvature_Component.h"
typedef boost::shared_ptr<Curvature_Component> Curvature_ComponentPtr;

double r_initial = 0.001 ;

template<typename T>
void writeTo(const std::string& filepath,const  vector<T>& data)
{

 ofstream filestream("C:\\vector.txt");
 std::copy(data.begin(),data.end(),std::ostream_iterator<T>( filestream," "));
 filestream.close();

}

void mepp_component_VSA_plugin::OnMouseLeftDown(QMouseEvent *event)
{
    cout<<"vertex selection done"<<endl;
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

//		if (doesExistComponentForViewer<VSA_ComponentPtr, VSA_Component>(viewer, polyhedron_ptr)) // important !!!
//		{
			mw->statusBar()->showMessage(tr("mepp_component_VSA_plugin: OnMouseLeftDown"), 1000);

			VSA_ComponentPtr component_ptr = findOrCreateComponentForViewer<VSA_ComponentPtr, VSA_Component>(viewer, polyhedron_ptr);

			if (viewer->getScenePtr()->get_loadType() == Time)
			{
				PolyhedronPtr new_polyhedron_ptr(new Polyhedron(*(viewer->getScenePtr()->get_polyhedron())));

				component_ptr->Get_Clicked_Vertices(viewer->getScenePtr()->get_polyhedron(), event->x(), event->y(), 10);
				viewer->recreateListsAndUpdateGL();
				SleeperThread::msleep(300);

				viewer->getScenePtr()->add_polyhedron(new_polyhedron_ptr);
				viewer->getScenePtr()->set_current_polyhedron(viewer->getScenePtr()->get_nb_polyhedrons()-1);

				viewer->setDynTitle();

				viewer->recreateListsAndUpdateGL();
			}
			else if (viewer->getScenePtr()->get_loadType() == Normal)
			{
				component_ptr->Get_Clicked_Vertices(viewer->getScenePtr()->get_polyhedron(), event->x(), event->y(), 10);
				viewer->recreateListsAndUpdateGL();
			}
//		}
	}
}

void mepp_component_VSA_plugin::OnKeyPress(QKeyEvent *event)
{
    double step_pct = 0.001;
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

    if ((event->key()==Qt::Key_H) & (double(r_initial+= step_pct) <=100.0))
        r_initial += step_pct ;
    if ((event->key()==Qt::Key_J) & (double(r_initial-= step_pct) >=0.1))
        r_initial -= step_pct ;

	if (doesExistComponentForViewer<VSA_ComponentPtr, VSA_Component>(viewer, polyhedron_ptr))
    {
    Viewer *viewer = (Viewer *)mw->activeMdiChild();// avoid bug under Linux
    PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
    VSA_ComponentPtr component_ptr = findOrCreateComponentForViewer<VSA_ComponentPtr, VSA_Component>(viewer, polyhedron_ptr);
    double maxDimInner_Val = component_ptr->getMaxDim(polyhedron_ptr);
    double r = (r_initial * maxDimInner_Val);
    component_ptr->draw_Sphere(polyhedron_ptr,(float)r, (int)1);
	viewer->recreateListsAndUpdateGL();
//				statusBar()->showMessage(tr("File saved: %1").arg(strippedName(fileName)), 2000);
    mw->statusBar()->showMessage(tr("Radius = %1 %").arg(r_initial*100.0), 1000);
	}
	}
}

void mepp_component_VSA_plugin::BurstWindSegmentation()
{
  if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer_org = (Viewer *)mw->activeMdiChild();


		emit(mw->get_actionNewEmpty()->trigger());

		QList<QMdiSubWindow *> lwindow = mw->mdiArea->subWindowList();

		for (int i=0; i<lwindow.size(); i++) // all viewers
		{
			Viewer* viewer = (Viewer *)qobject_cast<QWidget *>(lwindow[i]->widget());
			if (viewer->getScenePtr()->get_polyhedron()->empty())
			{
				PolyhedronPtr new_polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
				PolyhedronPtr polyhedron_ptr = viewer_org->getScenePtr()->get_polyhedron();

				QApplication::setOverrideCursor(Qt::WaitCursor);
				new_polyhedron_ptr->copy_from(&(*polyhedron_ptr));

				QApplication::restoreOverrideCursor();

				viewer->showAllScene();

				viewer->getScenePtr()->setcurrentFile(viewer_org->userFriendlyCurrentFile());
				viewer->setDynTitle();

				QApplication::setOverrideCursor(Qt::WaitCursor);

				/* Attention : tout changement d'un des paramètres ou des approches de notre algo exige refaire des tests pour ajuster les autres paramètres
				*/

				 int iter = 5;
                 int i = 0;
                 double sumCot = 0;
                 bool IsGeo= true;
                 double RadiusVal=0.01;


                         // we use CGAL_Example component here (as usual)
                VSA_ComponentPtr component_ptr = findOrCreateComponentForViewer<VSA_ComponentPtr, VSA_Component>(viewer, polyhedron_ptr);
                double maxdim=component_ptr->getMaxDim(polyhedron_ptr);
                Curvature_ComponentPtr component_ptr_curvature = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);
                 component_ptr->computeGaussianCurv1Ring(polyhedron_ptr);
                //component_ptr_curvature->principal_curvature(polyhedron_ptr,IsGeo,RadiusVal*maxdim);

                double sumCot2= MedianofVector(component_ptr->computeLaplacianCotg3(polyhedron_ptr));
                cout << "smoothing start" << endl;

                do
                {
                 sumCot = sumCot2;
                 polyhedron_ptr = component_ptr->AdatedVersionOfTaubinSmoothing(polyhedron_ptr,0.5, iter, true);
                 component_ptr->computeGaussianCurv1Ring(polyhedron_ptr);
                 //component_ptr_curvature->principal_curvature(polyhedron_ptr,IsGeo,RadiusVal*maxdim);
                 sumCot2= MedianofVector(component_ptr->computeLaplacianCotg3(polyhedron_ptr));
                 i= i+1;
                } while ((sumCot - sumCot2) > 0.0001);


               cout << "smoothing done" << endl;
               cout << "i =" <<i << endl;
                component_ptr->set_init(2); //(7)
                viewer->recreateListsAndUpdateGL();




              std::vector<double> dev = component_ptr->DeviationPolyhedronToPolyhedron(new_polyhedron_ptr, polyhedron_ptr);
              new_polyhedron_ptr = component_ptr->DepthFacets(new_polyhedron_ptr);
              dev.clear();

              for (Facet_iterator pFacet = new_polyhedron_ptr->facets_begin(); pFacet != new_polyhedron_ptr->facets_end(); pFacet++)
                {
                    dev.push_back(pFacet->depthfacet);

                }

            //  component_ptr->ConstructFaceDepthMap(new_polyhedron_ptr);

//            std::vector<double> grad= component_ptr->computeGradientMagnitudeFacet(new_polyhedron_ptr);
//                new_polyhedron_ptr = component_ptr->GradientDepthFacets(new_polyhedron_ptr);
//                cout << "done" << endl;
//                ConstructFaceGradientDepthMap(new_polyhedron_ptr);
//                //threshold =  MeanofVector(grad);
//
//                //threshold = component_ptr->getOtsuThreshold(grad);
//
           //  component_ptr->Burst_Wind_Segmentation_ON_Facets_2(new_polyhedron_ptr, dev,grad);


             //component_ptr->ConstructFaceColorMap(new_polyhedron_ptr);
//                 //new_polyhedron_ptr= component_ptr->MedianGradientSmoothing(new_polyhedron_ptr);
//
//                //ConstructGradientCotColorMap (new_polyhedron_ptr,threshold);
////                  cout << "thresold ::::::::::::::::" << threshold << endl;
////                WriteTxt2(grad, "/home/isit/cotang.txt");
////                component_ptr->set_init(3);
////                viewer->recreateListsAndUpdateGL();


			}
		}
	}
}

void mepp_component_VSA_plugin::BurstWindSegmentation2()
{
  if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer_org = (Viewer *)mw->activeMdiChild();


		emit(mw->get_actionNewEmpty()->trigger());

		QList<QMdiSubWindow *> lwindow = mw->mdiArea->subWindowList();

		for (int i=0; i<lwindow.size(); i++) // all viewers
		{
			Viewer* viewer = (Viewer *)qobject_cast<QWidget *>(lwindow[i]->widget());
			if (viewer->getScenePtr()->get_polyhedron()->empty())
			{
				PolyhedronPtr new_polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
				PolyhedronPtr polyhedron_ptr = viewer_org->getScenePtr()->get_polyhedron();

				QApplication::setOverrideCursor(Qt::WaitCursor);
				new_polyhedron_ptr->copy_from(&(*polyhedron_ptr));

				QApplication::restoreOverrideCursor();

				viewer->showAllScene();

				viewer->getScenePtr()->setcurrentFile(viewer_org->userFriendlyCurrentFile());
				viewer->setDynTitle();

				QApplication::setOverrideCursor(Qt::WaitCursor);

				/* Attention : tout changement d'un des paramètres ou des approches de notre algo exige refaire des tests pour ajuster les autres paramètres
				*/

				 int iter = 5;
                 int i = 0;
                 double sumCot = 0;
                 bool IsGeo= true;
                 double RadiusVal=r_initial;


                         // we use CGAL_Example component here (as usual)



                VSA_ComponentPtr component_ptr = findOrCreateComponentForViewer<VSA_ComponentPtr, VSA_Component>(viewer_org, polyhedron_ptr);
                double maxdim= component_ptr->getMaxDim(polyhedron_ptr);

                Curvature_ComponentPtr component_ptr_curvature = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer_org, polyhedron_ptr);
               //component_ptr->computeGaussianCurv1Ring(polyhedron_ptr);
                 component_ptr_curvature->principal_curvature(polyhedron_ptr,IsGeo,RadiusVal*maxdim);


                std::vector<double> curv= component_ptr->computeLaplacianCotg2(polyhedron_ptr);
                double sumCot2= MedianofVector(curv);
//                std::vector<double>::iterator it = std::max_element(curv.begin(), curv.end());
//                double sumCot2= *it;
                cout << "sumCot2 =" << sumCot2 <<endl;
                cout << "smoothing start" << endl;

//                do
//                {
                 sumCot = sumCot2;
                 polyhedron_ptr = component_ptr->Multi_scaleLaplacianSmoothing(polyhedron_ptr,0.5, iter, true, RadiusVal*maxdim);
                 //component_ptr->computeGaussianCurv1Ring(polyhedron_ptr);
                 component_ptr_curvature->principal_curvature(polyhedron_ptr,IsGeo,RadiusVal*maxdim);

                 curv= component_ptr->computeLaplacianCotg2(polyhedron_ptr);
                 sumCot2= MedianofVector(curv);
//                  std::vector<double>::iterator it = std::max_element(curv.begin(), curv.end());
//                 sumCot2= *it;
                 cout << "sumCot2 =" << sumCot2 <<endl;
                 i= i+1;
                 cout << "sumCot - sumCot2" << sumCot - sumCot2 << endl;
  //              } while (sumCot2 > 1e-8);
                // la condition de lissage qui fait le moins d'iteration de lissage et qui donnedes bonnes résultats 0.001

               cout << "smoothing done" << endl;
               cout << "i =" <<i << endl;
                component_ptr->set_init(2); //(7)
                viewer->recreateListsAndUpdateGL();
                viewer_org->recreateListsAndUpdateGL();




              std::vector<double> dev = component_ptr->DeviationPolyhedronToPolyhedron(new_polyhedron_ptr, polyhedron_ptr);
//              cout << "nous sommes la " << endl;
//              new_polyhedron_ptr = component_ptr->TranslateToZeroTheDepth(new_polyhedron_ptr,dev);
//              new_polyhedron_ptr = component_ptr->AdjustDepth(new_polyhedron_ptr);
              new_polyhedron_ptr = component_ptr->DepthFacets(new_polyhedron_ptr);
              dev.clear();

              for (Facet_iterator pFacet = new_polyhedron_ptr->facets_begin(); pFacet != new_polyhedron_ptr->facets_end(); pFacet++)
                {
                    dev.push_back(pFacet->depthfacet);

                }

             component_ptr->ConstructFaceDepthMap(new_polyhedron_ptr);
//
           std::vector<double> grad= component_ptr->computeGradientMagnitudeFacet(&new_polyhedron_ptr);
           double threshold= component_ptr->MedianofVector(grad);
           cout << "threshold =" << threshold << endl;
           component_ptr->ConstructGradFaceDepthMap(new_polyhedron_ptr, threshold);

//
           new_polyhedron_ptr= component_ptr->Burst_Wind_Segmentation_ON_Facets_2(new_polyhedron_ptr, dev, grad);
           cout << "nous sommes la " << endl;
//          // new_polyhedron_ptr = component_ptr->Morpho_Fill_holes(new_polyhedron_ptr);
//
//
//
//           //the minimum area of the segmented regions
           std::vector<double> areas;

          double minarea =component_ptr->m_Table_Proxy[0].Area;
          areas.push_back(minarea);
          double area;
          for(int i=1;i<component_ptr->m_Table_Proxy.size();i++)//For each proxy we select triangles that could grow
            {
                area= component_ptr->m_Table_Proxy[i].Area;
                areas.push_back(area);
                if (area < minarea)
                minarea= area;

            }


          double meanarea= mean_of(areas);
          sort(areas.begin(), areas.end());
          double epsilon = areas[areas.size()*1/4];
          double coeff;
//
           for(int i=0;i<component_ptr->m_Table_Proxy.size();i++)//For each proxy we select triangles that could grow
                {

                         area= component_ptr->m_Table_Proxy[i].Area;
                         RadiusVal= std::sqrt(area);

                         if (area < epsilon)
                         {
                            coeff = 1;

                         }
                         else if (area > epsilon && area < 2*epsilon)
                         {
                          coeff = 0.5;

                         }
                         else if (area > 2*epsilon && area < 3*epsilon)
                         {
                          coeff = 0.3;
                         }
                         else
                         {
                         coeff = 0.1;
                         }
                         component_ptr->closing(new_polyhedron_ptr,i, coeff*RadiusVal);
                        //new_polyhedron_ptr= component_ptr->Morpho_Fill_holes(new_polyhedron_ptr,i);
                         for (Facet_iterator pFacet = new_polyhedron_ptr->facets_begin(); pFacet != new_polyhedron_ptr->facets_end(); pFacet++)
                        {
                            if  (pFacet->morpho_closed !=0)
                                 pFacet->LabelVSA=i;
                        }

                 }


          write_DepthReg_in_Text("depth_armadillo.txt", new_polyhedron_ptr);
          write_LabelReg_in_Text("label_armadillo.txt", new_polyhedron_ptr);
//
            component_ptr->ConstructFaceColorMap(new_polyhedron_ptr);
//                 //new_polyhedron_ptr= component_ptr->MedianGradientSmoothing(new_polyhedron_ptr);
//
//                //ConstructGradientCotColorMap (new_polyhedron_ptr,threshold);
////                  cout << "thresold ::::::::::::::::" << threshold << endl;
////                WriteTxt2(grad, "/home/isit/cotang.txt");
////                component_ptr->set_init(3);
////                viewer->recreateListsAndUpdateGL();


			}
		}
	}
}

void mepp_component_VSA_plugin::write_DepthReg_in_Text(string output_name, PolyhedronPtr pMesh)
{
			std::ofstream file(output_name.c_str());
			for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
			{
                file << pFacet->depthfacet << endl;
			}
			file.flush();
			file.close();
}


void mepp_component_VSA_plugin::write_LabelReg_in_Text(string output_name, PolyhedronPtr pMesh)
{
			std::ofstream file(output_name.c_str());
			for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
			{
                file << pFacet->LabelVSA << endl;
			}
			file.flush();
			file.close();
}

double mepp_component_VSA_plugin::mean_of(const std::vector<double>& values)
{
  double sum = 0.0;
  for(unsigned long int i = 0; i < values.size(); ++i)
    sum += values[i];
  return sum / static_cast<double>(values.size());
}

double mepp_component_VSA_plugin::variance_of(const std::vector<double>& values)
{
  double mean = mean_of(values);
  double sum = 0.0;
  for(unsigned long int i = 0; i < values.size(); ++i)
    sum += (values[i] - mean) * (values[i] - mean);
  return sum / static_cast<double>(values.size() - 1);
}

// Fonction calculant l'écart-type des nombres contenus dans un vecteur
double mepp_component_VSA_plugin::standard_deviation_of(const std::vector<double>& values)
{
  return sqrt(variance_of(values));
}

double mepp_component_VSA_plugin::MedianofVector (std::vector<double> grades)
{

    sort(grades.begin(), grades.end());
    double median= *(grades.begin()+grades.size()/2); //89

    return median;
}

void mepp_component_VSA_plugin::FaceLabelsToColorMap()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		VSA_ComponentPtr component_ptr = findOrCreateComponentForViewer<VSA_ComponentPtr, VSA_Component>(viewer, polyhedron_ptr);
		component_ptr->ConstructFaceColorMap(polyhedron_ptr);

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_VSA_plugin, mepp_component_VSA_plugin);
#endif

#endif
