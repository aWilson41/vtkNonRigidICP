#include "vtkNonRigidICP.h"
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSTLReader.h>

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "   Error, Usage: vtkNonRigidICP <source stl file> <target stl file>" << std::endl;
		return EXIT_FAILURE;
	}

	// Input parameters
	std::string targetFile(argv[1]);
	std::string sourceFile(argv[2]);


	// Read in input data
	vtkNew<vtkSTLReader> targetReader;
	targetReader->SetFileName(targetFile.c_str());
	targetReader->Update();
	vtkNew<vtkSTLReader> sourceReader;
	sourceReader->SetFileName(sourceFile.c_str());
	sourceReader->Update();


	// Compute non rigid icp
	vtkNew<vtkNonRigidICP> surfaceRegister;
	surfaceRegister->SetMovingPolyData(sourceReader->GetOutput());
	surfaceRegister->SetFixedPolyData(targetReader->GetOutput());
	surfaceRegister->SetStiffness(1000.0);
	surfaceRegister->SetNumberOfIterations(40);
	surfaceRegister->SetOutputDisplacements(true);
	surfaceRegister->Update();
	vtkPolyData* outputPoly = surfaceRegister->GetOutput();


	// Setup render
	vtkDataArray* scalarArray = outputPoly->GetPointData()->GetScalars();
	scalarArray->SetName("Scalars");
	outputPoly->GetPointData()->SetActiveVectors("Scalars");

	// Setup render items
	vtkNew<vtkArrowSource> glyphSource;
	glyphSource->SetShaftRadius(0.05);
	glyphSource->SetTipRadius(0.15);
	glyphSource->SetTipLength(0.35);
	glyphSource->Update();
	vtkNew<vtkGlyph3D> glyphFilter;
	glyphFilter->SetInputData(outputPoly);
	glyphFilter->SetSourceData(glyphSource->GetOutput());
	glyphFilter->OrientOn();
	glyphFilter->SetScaleFactor(1.0);
	glyphFilter->SetScaleModeToScaleByVector();
	glyphFilter->SetColorModeToColorByScale();
	glyphFilter->Update();
	vtkNew<vtkPolyDataMapper> glyphMapper;
	glyphMapper->SetInputData(glyphFilter->GetOutput());
	glyphMapper->Update();
	vtkNew<vtkActor> displacementArrowsActor;
	displacementArrowsActor->SetMapper(glyphMapper);

	vtkNew<vtkPolyDataMapper> fixedPolyMapper;
	fixedPolyMapper->SetInputData(targetReader->GetOutput());
	fixedPolyMapper->Update();
	vtkNew<vtkActor> fixedPolyActor;
	fixedPolyActor->SetMapper(fixedPolyMapper);

	vtkNew<vtkPolyDataMapper> movingPolyMapper;
	movingPolyMapper->SetInputData(sourceReader->GetOutput());
	movingPolyMapper->Update();
	vtkNew<vtkActor> movingPolyActor;
	movingPolyActor->SetMapper(movingPolyMapper);
	movingPolyActor->GetProperty()->SetDiffuseColor(0.8, 0.2, 0.2);
	movingPolyActor->GetProperty()->SetOpacity(0.5);


	vtkNew<vtkRenderWindow> renWin;
	vtkNew<vtkRenderer> ren;
	renWin->AddRenderer(ren);
	vtkNew<vtkRenderWindowInteractor> iren;
	iren->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
	iren->SetRenderWindow(renWin);

	ren->SetBackground(0.75, 0.55, 0.62);
	ren->AddActor(displacementArrowsActor);
	ren->AddActor(fixedPolyActor);
	ren->AddActor(movingPolyActor);

	renWin->Render();
	iren->Start();

	return EXIT_SUCCESS;
}
