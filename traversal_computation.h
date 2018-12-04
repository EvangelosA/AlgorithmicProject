#ifndef TRAVERSAL_COMPUTATION__H__
#define TRAVERSAL_COMPUTATION__H__

typedef struct indx
{
	int i;
	int j;
}indx;

typedef struct traversalArrayInfo
{
	indx *traversalArray;
	int traversalArraySize;
	int nextAvailablePos;
}traversalArrayInfo;

traversalArrayInfo* DFD_traversal(curveInfo, curveInfo);
curveInfo *meanCurve(curveInfo *, curveInfo *);

#endif