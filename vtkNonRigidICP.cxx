#include "vtkNonRigidICP.h"
#include <Eigen/Sparse>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
//#include <chrono>

vtkStandardNewMacro(vtkNonRigidICP);

using Triplet = Eigen::Triplet<double>;

void vtkNonRigidICP::Cleanup()
{
	outputCellArray = NULL;
	fixedKdTree = NULL;
}

int vtkNonRigidICP::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVec, vtkInformationVector* outputVec)
{
	// Get input polys
	vtkInformation* inInfo1 = inputVec[0]->GetInformationObject(0);
	vtkPolyData* fixedPolyData = vtkPolyData::SafeDownCast(inInfo1->Get(vtkDataObject::DATA_OBJECT()));
	vtkInformation* inInfo2 = inputVec[1]->GetInformationObject(0);
	vtkPolyData* movingPolyData = vtkPolyData::SafeDownCast(inInfo2->Get(vtkDataObject::DATA_OBJECT()));
	// Get output poly
	vtkInformation* outInfo = outputVec->GetInformationObject(0);
	vtkPolyData* outputMovingPolyData = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	outputMovingPolyData->DeepCopy(movingPolyData);

	// Only compute if we haven't stored a previous one
	if (outputCellArray == NULL)
	{
		// Compute the edge only topology of the output polygon
		outputCellArray = vtkSmartPointer<vtkCellArray>::New();
		for (vtkIdType i = 0; i < outputMovingPolyData->GetNumberOfCells(); i++)
		{
			vtkCell* cell = outputMovingPolyData->GetCell(i);
			if (cell->GetNumberOfPoints() == 3)
			{
				vtkIdType ptIds1[2] = { cell->GetPointId(0), cell->GetPointId(1) };
				outputCellArray->InsertNextCell(2, ptIds1);
				vtkIdType ptIds2[2] = { cell->GetPointId(1), cell->GetPointId(2) };
				outputCellArray->InsertNextCell(2, ptIds2);
				vtkIdType ptIds3[2] = { cell->GetPointId(2), cell->GetPointId(0) };
				outputCellArray->InsertNextCell(2, ptIds3);
			}
			else if (cell->GetNumberOfPoints() == 2)
			{
				vtkIdType ptIds[2] = { cell->GetPointId(0), cell->GetPointId(1) };
				outputCellArray->InsertNextCell(2, ptIds);
			}
		}
	}
	// Setup a kdtree for nearest neighbor searches on the fixed poly
	if (fixedKdTree == NULL)
	{
		fixedKdTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
		fixedKdTree->SetDataSet(fixedPolyData);
		fixedKdTree->BuildLocator();
	}

	double* targetBounds = fixedPolyData->GetBounds();
	double maxDist = std::pow(std::pow(targetBounds[1] - targetBounds[0], 2) + std::pow(targetBounds[3] - targetBounds[2], 2) + std::pow(targetBounds[5] - targetBounds[4], 2), 1.0 / 3.0);

	const vtkIdType fixedPolyNumberOfPoints = fixedPolyData->GetNumberOfPoints();
	const vtkIdType movingPolyNumberOfPoints = outputMovingPolyData->GetNumberOfPoints();
	const vtkIdType movingPolyNumberOfCells = outputCellArray->GetNumberOfCells();

	double stiffness = Stiffness;
	vtkNew<vtkIdTypeArray> correspondences;
	correspondences->SetName("Correspondences");
	correspondences->SetNumberOfValues(movingPolyNumberOfPoints);

	for (unsigned int iter = 0; iter < NumberOfIterations; iter++)
	{
		vtkNew<vtkKdTreePointLocator> movingKdTree;
		movingKdTree->SetDataSet(movingPolyData);
		movingKdTree->Update();

		vtkNew<vtkPolyDataNormals> normalsFilter;
		normalsFilter->SetInputData(outputMovingPolyData);
		normalsFilter->Update();
		vtkDataArray* normalsArray = normalsFilter->GetOutput()->GetPointData()->GetNormals();

		correspondences->Fill(-1);
		// Find the forward correspondences of all the moving points
		for (vtkIdType i = 0; i < movingPolyNumberOfPoints; i++)
		{
			double pt1[3];
			outputMovingPolyData->GetPoint(i, pt1);

			const vtkIdType closestPtIndex = fixedKdTree->FindClosestPoint(pt1);
			correspondences->SetValue(i, closestPtIndex);
		}

		// Construct matrix A, composed of MG and WD
		Eigen::SparseMatrix<double> A(4 * movingPolyNumberOfCells + movingPolyNumberOfPoints, 4 * movingPolyNumberOfPoints);

		std::vector<Triplet> ATripletList;
		ATripletList.reserve(movingPolyNumberOfCells * 8 + 4 * movingPolyNumberOfPoints);

		// Construct M kronecker G, regularization term resulting from the topology of the mesh
		vtkNew<vtkIdList> idList;
		outputCellArray->InitTraversal();
		for (vtkIdType i = 0; i < movingPolyNumberOfCells; i++)
		{
			outputCellArray->GetNextCell(idList);
			const vtkIdType id1 = idList->GetId(0);
			const vtkIdType id2 = idList->GetId(1);

			for (unsigned int j = 0; j < 3; j++)
			{
				ATripletList.push_back(Triplet(4 * i + j, 4 * id1 + j, 1.0  * stiffness));
			}
			ATripletList.push_back(Triplet(4 * i + 3, 4 * id1 + 3, SkewWeight * stiffness));

			for (unsigned int j = 0; j < 3; j++)
			{
				ATripletList.push_back(Triplet(4 * i + j, 4 * id2 + j, -1.0 * stiffness));
			}
			ATripletList.push_back(Triplet(4 * i + 3, 4 * id2 + 3, -SkewWeight * stiffness));
		}
		// Construct WD, for the distance term
		for (vtkIdType i = 0; i < movingPolyNumberOfPoints; i++)
		{
			const double w = (correspondences->GetValue(i) != -1) ? 1.0 : 0.0;
			double pt[3];
			outputMovingPolyData->GetPoint(i, pt);
			for (vtkIdType j = 0; j < 3; j++)
			{
				ATripletList.push_back(Triplet(4 * movingPolyNumberOfCells + i, i * 4 + j, pt[j] * w));
			}
			ATripletList.push_back(Triplet(4 * movingPolyNumberOfCells + i, i * 4 + 3, 1.0 * w));
		}
		A.setFromTriplets(ATripletList.begin(), ATripletList.end());

		Eigen::SparseMatrix<double> AT(A.cols(), A.rows());
		AT = A.transpose();
		Eigen::SparseMatrix<double> ATA(4 * movingPolyNumberOfPoints, 4 * movingPolyNumberOfPoints);
		ATA = AT * A;


		// Construct matrix B, resulting from correspodences of the points (data term)
		Eigen::MatrixXd B(4 * movingPolyNumberOfCells + movingPolyNumberOfPoints, 3);
		B.setZero();

		for (unsigned int i = 0; i < movingPolyNumberOfPoints; i++)
		{
			const vtkIdType correspondence = correspondences->GetValue(i);
			const double w = (correspondence != -1) ? 1.0 : 0.0;
			double pt[3];
			fixedPolyData->GetPoint(correspondence, pt);
			for (int j = 0; j < 3; j++)
			{
				B(4 * movingPolyNumberOfCells + i, j) = pt[j] * w;
			}
		}

		Eigen::MatrixXd ATB(4 * movingPolyNumberOfPoints, 3);
		ATB.setZero();
		ATB = AT * B;


		// Write the matrices
		/*Eigen::saveMarket(AmatrixFinal, "C:/Users/Andx_/Desktop/a.txt");
		Eigen::saveMarket(BmatrixFinal, "C:/Users/Andx_/Desktop/b.txt");*/

		Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
		solver.analyzePattern(ATA);
		solver.factorize(ATA);
		if (solver.info() != Eigen::Success)
			printf("Decomposition Failed.\n");
		Eigen::MatrixXd X(4 * movingPolyNumberOfPoints, 3);
		X.setZero();
		X = solver.solve(ATB);

		// Do the transformation
		vtkPoints* pts = outputMovingPolyData->GetPoints();
		for (vtkIdType i = 0; i < movingPolyNumberOfPoints; i++)
		{
			double pt[3];
			pts->GetPoint(i, pt);

			double results[3] = { 0.0, 0.0, 0.0 };
			for (int j = 0; j < 3; j++)
			{
				results[j] =
					pt[0] * X(i * 4 + 0, j) +
					pt[1] * X(i * 4 + 1, j) +
					pt[2] * X(i * 4 + 2, j) +
					1.0 * X(i * 4 + 3, j);
			}
			//printf("%f,%f,%f,%f\n", X(i * 4 + 0, 0), X(i * 4 + 1, 1), X(i * 4 + 2, 2), X(i * 4 + 3, 3));

			pts->SetPoint(i, results[0], results[1], results[2]);
		}

		// Invokes the progress update event
		if (iter % 4 == 0)
			UpdateProgress(static_cast<double>(iter) / NumberOfIterations);
	}

	if (!InteractiveExecute)
		Cleanup();

	if (OutputDisplacements)
	{
		// Stash a copy of the registered/moved surface
		vtkNew<vtkPolyData> registeredPolyData;
		registeredPolyData->DeepCopy(outputMovingPolyData);

		// Bring back the original polygon
		outputMovingPolyData->DeepCopy(movingPolyData);

		// Setup an array of displacements
		vtkNew<vtkDoubleArray> displacements;
		displacements->SetName("Displacements");
		displacements->SetNumberOfComponents(3);
		displacements->SetNumberOfTuples(movingPolyNumberOfPoints);

		// The indices should be the same so just compute vector differences
		for (vtkIdType i = 0; i < movingPolyNumberOfPoints; i++)
		{
			double orgPt[3];
			double movedPt[3];
			registeredPolyData->GetPoint(i, movedPt);
			outputMovingPolyData->GetPoint(i, orgPt);
			displacements->SetTuple3(i, movedPt[0] - orgPt[0], movedPt[1] - orgPt[1], movedPt[2] - orgPt[2]);
		}

		outputMovingPolyData->GetPointData()->SetScalars(displacements);
	}

	return 1;
}