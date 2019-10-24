#pragma once
#include <vtkPolyDataAlgorithm.h>

class vtkCellArray;
class vtkKdTreePointLocator;

/// Description
// This filter registers two vtkPolyData's iteratively. Can either run in one go or from an old
// state for later runs. Useful for computing one iteration, pausing, then computing the next later...
// However, it requires extra data (a spatial partitioning) to be stored for efficieny which can be
// cleaned with Cleanup()
/// Parameters
// The stiffness regularizes the transformation by smoothening transforms with neighboring transforms
// The skew weight controls the type of transformations allowed. A larger value will favor skews letting
// the points more easily slide past each other.
class vtkNonRigidICP : public vtkPolyDataAlgorithm
{
public:
	static vtkNonRigidICP* New();
	vtkTypeMacro(vtkNonRigidICP, vtkPolyDataAlgorithm);

	vtkNonRigidICP() { SetNumberOfInputPorts(2); }

	void SetFixedPolyData(vtkPolyData* data) { SetInputData(0, data); }
	void SetMovingPolyData(vtkPolyData* data) { SetInputData(1, data); }

	vtkGetMacro(Epsilon, double);
	vtkSetMacro(Epsilon, double);

	vtkGetMacro(NumberOfIterations, unsigned int);
	vtkSetMacro(NumberOfIterations, unsigned int);

	vtkGetMacro(Stiffness, double);
	vtkSetMacro(Stiffness, double);

	vtkGetMacro(SkewWeight, double);
	vtkSetMacro(SkewWeight, double);

	vtkGetMacro(InteractiveExecute, bool);
	vtkSetMacro(InteractiveExecute, bool);

	vtkGetMacro(OutputDisplacements, bool);
	vtkSetMacro(OutputDisplacements, bool);

	// If InteractiveExecute is off this is called at the end of the update
	// otherwise it is called either on change of inputs or user can call.
	void Cleanup();

protected:
	int RequestData(vtkInformation* request, vtkInformationVector** inputVec, vtkInformationVector* outputVec) VTK_OVERRIDE;

private:
	vtkNonRigidICP(const vtkNonRigidICP&);
	void operator=(const vtkNonRigidICP&);

private:
	vtkSmartPointer<vtkCellArray> outputCellArray;
	vtkSmartPointer<vtkKdTreePointLocator> fixedKdTree;

private:
	unsigned int NumberOfIterations = 20;
	double Epsilon = 0.1;
	double Stiffness = 1000.0;
	double SkewWeight = 1.0;
	bool OutputDisplacements = false;
	bool InteractiveExecute = false;
};