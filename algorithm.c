#include<stdio.h>
#include<windows.h>
#include<stdlib.h>
#include <limits.h>
#include"stackL.h"
#include"adjList.h"
#include"DFS.h"
#include"LinkedQueue.h"
#include"BFS.h"
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define ROW 5
#define COL 5
#define TRUE 1
#define FALSE 0
#define MAX_VERTICES 100
#define MAX_EDGES 1000
#define MAX 100
#define MAX_QUEUE_SIZE 10
#define INF 1000 // 무한대
#define MAX_N 10

int inputArray[MAX_N];   // 입력 배열
int result[MAX_N]; // 조합을 저장할 배열
int isVisited[MAX_N]; // 방문 여부를 체크할 배열
int sorted[MAX];
extern int size;
int temp_arr[MAX];
int graph[MAX_VERTICES][MAX_EDGES];
int visited[MAX_VERTICES];
void selectSort(int arr[], int size) {

	int i, j, t, min, temp;

	for (i = 0; i < size - 1; i++) {
		min = i;
		for (j = i + 1; j < size; j++) {
			if (arr[j] < arr[min]) min = j;
		}

		temp = arr[i];
		arr[i] = arr[min];
		arr[min] = temp;
		printf("\n%d단계 : ", i + 1);
		for (t = 0; t < size; t++) {
			printf("%3d", arr[t]);
		}
	}
}
void BubbleSort(int arr[], int size) {
	int i, j, t, temp;

	for (i = size - 1; i > 0; i--) {
		printf("\n %d단계>> ", size - i);
		for (j = 0; j < i; j++) {
			if (arr[j] > arr[j + 1]) {
				temp = arr[j];
				arr[j] = arr[j + 1];
				arr[j + 1] = temp;
			}

			printf("\n\t");
			for (t = 0; t < size; t++) {
				printf("%3d ", arr[t]);
			}
		}
	}
}
void InsertionSort(int arr[], int size) {
	int i, j, t, temp;

	for (i = 1; i < size; i++) {
		temp = arr[i];
		j = i;
		while ((j > 0) && (arr[j - 1] > temp)) {
			arr[j] = arr[j - 1];
			j = j - 1;
		}
		arr[j] = temp;
		printf("\n %d단계 : ", i);
		for (t = 0; t < size; t++) printf("%3d ", arr[t]);
	}
}
void merge(int* arr, int start, int middle, int end) {
	int t;
	int i = start;
	int j = middle + 1;
	int k = start;

	//비교하여 데이터정렬 및 삽입
	while (i <= middle && j <= end) {
		if (arr[i] <= arr[j]) sorted[k++] = arr[i++];
		else  sorted[k++] = arr[j++];
	
	}

	//남은 데이터 삽입
	if (i > middle) {
		for (t = j; t <= end; t++, k++) {
			sorted[k] = arr[t];
		}
	}
	else {
		for (t = i; t <= middle; t++, k++) {
			sorted[k] = arr[t];
		}
	}

	for (t = start; t <= end; t++) {
		arr[t] = sorted[t];
	}

	printf("\n 병합 정렬 >> ");
	for (int t = 0; t <= end; t++) {
		printf("%4d ", arr[t]);
	}
}
void mergeSort(int* arr, int start, int end) {
	//크기가 1 일대 까지 호출, 1단위 까지 쪼갬
	if (start < end) {
		int middle = (start + end) / 2;
		mergeSort(arr, start, middle);
		mergeSort(arr, middle + 1, end);
		//다시 병합
		merge(arr, start, middle, end);
	}
}
int i = 0;
int partition(int arr[], int begin, int end, int size) {
	int pivot, L, R, t, temp;
	L = begin;
	R = end;
	pivot = (begin + end) / 2;
	printf("\n [%d단계 : pivot = %d ] \n", ++i, arr[pivot]);
	while (L < R) {
		while ((arr[L] < arr[pivot]) && (L < R)) L++;
		while ((arr[R] >= arr[pivot]) && (L < R)) R--;
		if (L < R) {
			temp = arr[L];
			arr[L] = arr[R];
			arr[R] = temp;

			if (L == pivot) pivot = R;
		}

	}

	temp = arr[pivot];
	arr[pivot] = arr[R];
	arr[R] = temp;
	for (t = 0; t < size; t++) printf("%4d ", arr[t]);
	printf("\n");
	return R;
}
void quickSort(int arr[], int begin, int end, int size) {
	int p;
	if (begin < end) {
		p = partition(arr, begin, end, size);
		quickSort(arr, begin, p - 1, size);
		quickSort(arr, p + 1, end, size);
	}
}
void printArray(int arr[], int size) {
	for (int t = 0; t < size; t++) {
		printf("%4d ", arr[t]);
	}
	printf("\n");
}
void HeapSort(int arr[], int size) {
	// Build Max Heap
	for (int i = 1; i < size; i++) {
		int c = i;
		do {
			int parent = (c - 1) / 2;
			if (arr[parent] < arr[c]) {
				int temp = arr[parent];
				arr[parent] = arr[c];
				arr[c] = temp;
			}
			c = parent;
		} while (c != 0);

		printf("빌드 후 최대 힙 ( %d): ", i);
		printArray(arr, size);
	}

	// Heap Sort
	for (int i = size - 1; i >= 0; i--) {
		int temp = arr[0];
		arr[0] = arr[i];
		arr[i] = temp;
		int parent = 0;
		int c = 1;

		do {
			c = parent * 2 + 1;

			if (arr[c] < arr[c + 1] && c < i - 1) {
				c++;
			}

			if (arr[parent] < arr[c] && c < i) {
				temp = arr[parent];
				arr[parent] = arr[c];
				arr[c] = temp;
			}

			parent = c;
		} while (c < i);

		printf("힙 지정 후 (단계 %d): ", size - i);
		printArray(arr, size);
	}
}
void radixSort(int* arr, int size) {
	int result[MAX], maxValue = 0;
	int exp = 1;

	// Find the maximum value in the array
	for (int i = 0; i < size; i++) {
		if (arr[i] > maxValue) maxValue = arr[i];
	}

	while (maxValue / exp > 0) { // Iterate through each digit place
		int bucket[10] = { 0 };

		// Counting the occurrences of each digit
		for (int i = 0; i < size; i++) {
			bucket[arr[i] / exp % 10]++;
		}

		// Accumulate the counts to determine the positions
		for (int i = 1; i < 10; i++) {
			bucket[i] += bucket[i - 1];
		}

		// Place the elements in sorted order based on the current digit
		for (int i = size - 1; i >= 0; i--) {
			result[--bucket[arr[i] / exp % 10]] = arr[i];
		}

		// Update the original array with the sorted order
		for (int i = 0; i < size; i++) {
			arr[i] = result[i];
		}

		// Print the current state of the array after each pass

		printf("%d단계 : ", 1 + i++);
		printArray(arr, size);
		exp *= 10;
	}
}
void countingSort(int* arr, int size) {
	int count[MAX];
	for (int i = 0; i < MAX; i++) {
		count[i] = 0; // 초기화 
	}

	for (int i = 0; i < size; i++) {
		int val = arr[i];
		count[val]++; // 개수 세기 
	}

	// Print the sorted array after each pass
	for (int i = 0; i < MAX; i++) {
		for (int j = 0; j < count[i]; j++) {
			printf("%d ", i);
		}

		// Print the current state of the array after each pass
		if (count[i] > 0) {
			printf(" (정렬 후 %d): ", i);
			for (int k = 0; k < size; k++) {
				if (arr[k] == i) {
					printf("%d ", i);
				}
			}
			printf("\n");
		}
	}
}
void initializeGraph() {
	int i, j;
	for (i = 0; i < MAX_VERTICES; i++) {
		visited[i] = 0;
		for (j = 0; j < MAX_EDGES; j++) {
			graph[i][j] = 0;
		}
	}
}
int isStackEmpty(void) {
	if (top == NULL) return 1;
	else return 0;
}
void push(element item) {
	stackNode* temp = (stackNode*)malloc(sizeof(stackNode));

	temp->data = item;
	temp->link = top;
	top = temp;
}
element pop(void) {
	element item;
	stackNode* temp = top;

	if (isStackEmpty()) {
		printf("\n\n Stack is empty! \n");
		return 0;
	}
	else {
		item = temp->data;
		top = temp->link;
		free(temp);
		return item;
	}
}
element peek(void) {
	if (isStackEmpty()) {
		printf("\n\n Stack is empty! \n");
		return 0;
	}
	else {
		return(top->data);
	}
}
void printStack(void) {
	stackNode* p = top;
	printf("\n STACK [ ");
	while (p) {
		printf("%d ", p->data);
		p = p->link;
	}
	printf("] ");
}

void createGraph(graphType* g) {
	int v;
	g->n = 0;
	for (v = 0; v < MAX_VERTEX; v++) {
		g->adjList_H[v] = NULL;
		g->visited[v] = FALSE;
	}
}
void insertVertex(graphType* g, int v) {
	if (((g->n) + 1) > MAX_VERTEX) {
		printf("\n 그래프 정점의 개수를 초과하였습니다!");
		return;
	}
	g->n++;
}
void insertEdge(graphType* g, int u, int v) {
	graphNode* node;

	if (u >= g->n || v >= g->n) {
		printf("\n 그래프에 없는 정점입니다!");
		return;
	}
	node = (graphNode*)malloc(sizeof(graphNode));
	node->vertex = v;
	node->link = g->adjList_H[u];
	g->adjList_H[u] = node;
}
void print_adjList(graphType* g) {
	int i;
	graphNode* p;
	for (i = 0; i < g->n; i++) {
		printf("\n\t\t정점 %c의 인접 리스트", i + 65);
		p = g->adjList_H[i];
		while (p) {
			printf(" -> %c", p->vertex + 65);
			p = p->link;
		}
	}
}
void DFS_adjList(graphType* g, int v) {
	graphNode* w;
	top = NULL;
	push(v);
	g->visited[v] = TRUE;
	printf(" %c", v + 65);

	while (!isStackEmpty()) {
		w = g->adjList_H[v];
		while (w) {
			if (!g->visited[w->vertex]) {
				push(w->vertex);
				g->visited[w->vertex] = TRUE;
				printf(" %c", w->vertex + 65);
				v = w->vertex;
				w = g->adjList_H[v];
			}
			else w = w->link;
		}
		v = pop();
	}
}

LQueueType* createLinkedQueue(void) {
	LQueueType* LQ;
	LQ = (LQueueType*)malloc(sizeof(LQueueType));
	LQ->front = NULL;
	LQ->rear = NULL;
	return LQ;
}
int isLQEmpty(LQueueType* LQ) {
	if (LQ->front == NULL) {
		//printf(" Linked Queue is empty! ");
		return 1;
	}
	else return 0;
}
void enLQueue(LQueueType* LQ, element item) {
	QNode* newNode = (QNode*)malloc(sizeof(QNode));
	newNode->data = item;
	newNode->link = NULL;
	if (LQ->front == NULL) {
		LQ->front = newNode;
		LQ->rear = newNode;
	}
	else {
		LQ->rear->link = newNode;
		LQ->rear = newNode;
	}
}
element deLQueue(LQueueType* LQ) {
	QNode* old = LQ->front;
	element item;
	if (isLQEmpty(LQ)) return;
	else {
		item = old->data;
		LQ->front = LQ->front->link;
		if (LQ->front == NULL)
			LQ->rear = NULL;
		free(old);
		return item;
	}
}
element peekLQ(LQueueType* LQ) {
	element item;
	if (isLQEmpty(LQ)) return;
	else {
		item = LQ->front->data;
		return item;
	}
}
void printLQ(LQueueType* LQ) {
	QNode* temp = LQ->front;
	printf(" LInked Queue : [");
	while (temp) {
		printf("%3d", temp->data);
		temp = temp->link;
	}
	printf(" ]");
}

void BFS_adjList(graphType* g, int v) {
	graphNode* w;
	LQueueType* Q;
	Q = createLinkedQueue();
	g->visited[v] = TRUE;
	printf(" %c", v + 65);
	enLQueue(Q, v);

	while (!isLQEmpty(Q)) {
		v = deLQueue(Q);
		for(w = g->adjList_H[v]; w; w = w->link)
			if (!g->visited[w->vertex]) {
				g->visited[w->vertex] = TRUE;
				printf(" %c", w->vertex + 65);
				enLQueue(Q, w->vertex);
			}
	}

}

typedef struct {
	int w;
	int p;
} Eedge;

void prim(int G[MAX_VERTICES][MAX_VERTICES], int V, int r) {
	int S[MAX_VERTICES]; // 정점 집합
	Eedge d[MAX_VERTICES]; // 거리 및 부모 정보

	for (int u = 0; u < V; u++) {
		S[u] = 0; // S를 빈 집합으로 초기화
		d[u].w = INT_MAX;
		d[u].p = -1;
	}

	d[r].w = 0;

	while (1) {
		// 가장 작은 거리 값을 갖는 정점을 찾음
		int u = -1;
		for (int i = 0; i < V; i++) {
			if (!S[i] && (u == -1 || d[i].w < d[u].w)) {
				u = i;
			}
		}

		if (u == -1) {
			break; // 모든 정점이 MST에 포함됨
		}

		S[u] = 1; // 선택된 정점을 집합에 추가

		// 인접 정점의 거리 값을 업데이트
		for (int v = 0; v < V; v++) {
			if (!S[v] && G[u][v] < d[v].w) {
				d[v].w = G[u][v];
				d[v].p = u;
			}
		}
	}

	// MST의 간선을 출력
	printf("최소 신장 트리의 간선:\n");
	for (int i = 0; i < V; i++) {
		if (d[i].p != -1) {
			printf("(%d, %d) - 가중치: %d\n", d[i].p, i, d[i].w);
		}
	}
}

typedef struct {
	int start;
	int end;
	int weight;
} Edge;
typedef struct {
	int parent;
	int rank;
} Subset;
int find(Subset subsets[], int i) {
	if (subsets[i].parent != i) {
		subsets[i].parent = find(subsets, subsets[i].parent);
	}
	return subsets[i].parent;
}
void Union(Subset subsets[], int x, int y) {
	int rootX = find(subsets, x);
	int rootY = find(subsets, y);

	if (subsets[rootX].rank < subsets[rootY].rank) {
		subsets[rootX].parent = rootY;
	}
	else if (subsets[rootX].rank > subsets[rootY].rank) {
		subsets[rootY].parent = rootX;
	}
	else {
		subsets[rootY].parent = rootX;
		subsets[rootX].rank++;
	}
}
int compareEdges(const void* a, const void* b) {
	return ((Edge*)a)->weight - ((Edge*)b)->weight;
}
void kruskal(Edge edges[], int V, int E) {
	Edge result[MAX_VERTICES];
	Subset subsets[MAX_VERTICES];

	for (int v = 0; v < V; v++) {
		subsets[v].parent = v;
		subsets[v].rank = 0;
	}

	qsort(edges, E, sizeof(Edge), compareEdges);

	int t = 0; // 인덱스 변수 및 결과 트리에 포함된 간선 수

	for (int i = 0; t < V - 1 && i < E; i++) {
		Edge nextEdge = edges[i];

		int x = find(subsets, nextEdge.start);
		int y = find(subsets, nextEdge.end);

		if (x != y) {
			result[t++] = nextEdge;
			Union(subsets, x, y);
		}
	}

	// 결과 출력
	printf("최소 신장 트리의 간선:\n");
	for (int i = 0; i < t; i++) {
		printf("(%d, %d) - 가중치: %d\n", result[i].start, result[i].end, result[i].weight);
	}
}



typedef struct {
	int row, col;
} Point;

typedef struct {
	int f, g, h;
	Point parent;
} HeapNode;

void swap(HeapNode* x, HeapNode* y) {
	HeapNode temp = *x;
	*x = *y;
	*y = temp;
}

void heapify(HeapNode heap[], int size, int i) {
	int smallest = i;
	int left = 2 * i + 1;
	int right = 2 * i + 2;

	if (left < size && heap[left].f < heap[smallest].f)
		smallest = left;

	if (right < size && heap[right].f < heap[smallest].f)
		smallest = right;

	if (smallest != i) {
		swap(&heap[i], &heap[smallest]);
		heapify(heap, size, smallest);
	}
}

void insertHeap(HeapNode heap[], int* size, HeapNode newNode) {
	(*size)++;
	int i = *size - 1;
	heap[i] = newNode;

	while (i > 0 && heap[(i - 1) / 2].f > heap[i].f) {
		swap(&heap[i], &heap[(i - 1) / 2]);
		i = (i - 1) / 2;
	}
}

HeapNode extractMin(HeapNode heap[], int* size) {
	if (*size == 0) {
		fprintf(stderr, "Heap is empty\n");
		exit(EXIT_FAILURE);
	}

	if (*size == 1) {
		(*size)--;
		return heap[0];
	}

	HeapNode root = heap[0];
	heap[0] = heap[(*size) - 1];
	(*size)--;
	heapify(heap, *size, 0);

	return root;
}

bool isValid(int row, int col) {
	return (row >= 0) && (row < ROW) && (col >= 0) && (col < COL);
}

bool isUnBlocked(int grid[ROW][COL], int row, int col) {
	return grid[row][col] == 1; // 1은 길, 0은 벽
}

bool isDestination(int row, int col, Point dest) {
	return (row == dest.row) && (col == dest.col);
}

int heuristicValue(int row, int col, Point dest) {
	return abs(row - dest.row) + abs(col - dest.col);
}

void tracePath(Point parent[][COL], Point dest) {
	printf("최단 경로: ");
	int row = dest.row;
	int col = dest.col;
	Point path[ROW * COL];
	int pathLength = 0;

	while (!(parent[row][col].row == -1 && parent[row][col].col == -1)) {
		path[pathLength].row = row;
		path[pathLength].col = col;
		pathLength++;

		Point temp = parent[row][col];
		row = temp.row;
		col = temp.col;
	}

	for (int i = pathLength - 1; i >= 0; i--)
		printf("-> (%d, %d) ", path[i].row, path[i].col);
}

void aStarSearch(int grid[ROW][COL], Point src, Point dest) {
	if (!isValid(src.row, src.col) || !isValid(dest.row, dest.col)) {
		fprintf(stderr, "유효하지 않은 시작 또는 목적지\n");
		exit(EXIT_FAILURE);
	}

	if (!isUnBlocked(grid, src.row, src.col) || !isUnBlocked(grid, dest.row, dest.col)) {
		fprintf(stderr, "시작 또는 목적지가 벽에 막혀있음\n");
		exit(EXIT_FAILURE);
	}

	if (isDestination(src.row, src.col, dest)) {
		fprintf(stderr, "이미 목적지에 도착함\n");
		exit(EXIT_FAILURE);
	}

	bool closedList[ROW][COL];
	memset(closedList, false, sizeof(closedList));

	Point parent[ROW][COL];
	memset(parent, -1, sizeof(parent));

	int row, col;
	HeapNode heap[ROW * COL];
	int heapSize = 0;

	HeapNode startNode = { 0, 0, 0, src };
	insertHeap(heap, &heapSize, startNode);

	while (heapSize > 0) {
		HeapNode currentNode = extractMin(heap, &heapSize);
		row = currentNode.parent.row;
		col = currentNode.parent.col;
		closedList[row][col] = true;

		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				int newRow = row + i;
				int newCol = col + j;

				if (isValid(newRow, newCol) && isUnBlocked(grid, newRow, newCol) && !closedList[newRow][newCol]) {
					int gNew = currentNode.g + 1;
					int hNew = heuristicValue(newRow, newCol, dest);
					int fNew = gNew + hNew;

					if (heapSize == ROW * COL) {
						fprintf(stderr, "힙이 가득 참\n");
						exit(EXIT_FAILURE);
					}

					if (parent[newRow][newCol].row == -1 || fNew < heap[0].f) {
						HeapNode newNode = { fNew, gNew, hNew, {newRow, newCol} };
						insertHeap(heap, &heapSize, newNode);
						parent[newRow][newCol] = currentNode.parent;
					}
				}
			}
		}
	}

	if (!closedList[dest.row][dest.col]) {
		printf("목적지에 도착할 수 없음\n");
	}
	else {
		tracePath(parent, dest);
	}
}

// n개 중에서 r개를 선택하는 조합을 찾는 함수
void combination(int n, int r, int index) {
	if (r == 0) {
		// 조합이 완성되면 결과 출력
		for (int i = 0; i < index; ++i) {
			printf("%d ", result[i]);
		}
		printf("\n");
		return;
	}

	for (int i = 0; i < n; ++i) {
		if (!isVisited[i]) {
			isVisited[i] = 1;         // 해당 원소를 선택했다고 표시
			result[index] = inputArray[i]; // 결과 배열에 저장
			combination(n, r - 1, index + 1); // 재귀 호출
			isVisited[i] = 0;         // 백트래킹
		}
	}
}

int sortnumber = 0;
int* arr;

int size = 0;

void setting() {
	printf("\n▶ 배열의 크기를 입력해주세요(최대 100까지): ");
	scanf_s("%d", &size);

	arr = (int*)malloc(sizeof(int) * size);

	printf("\n▶ 배열을 %d개 입력해주세요(띄어쓰기로 구분): ", size);
	for (int i = 0; i < size; i++) {
		scanf_s("%d", &arr[i]);
	}
}
void termination() {
	int number = 0;
	printf("\n\n▷ 알고리즘 실행을 완료하였습니다.\n");
	printf("=================================================\n");
	printf("l                  0.메인 메뉴                  l\n");
	printf("=================================================\n");
	printf("l                  1.종료                       l\n");
	printf("=================================================\n\n");

	printf("▶ 번호를 입력하세요: ");
	scanf_s("%d", &number);

	if (number == 0) {
		system("cls");
		main();
	}
	else {
		printf("▷ 프로그램을 종료합니다.\n");
		return 0;
	}
}
int main() {
	while (1) {
		int number;
		printf("\n");
		printf("\n");
		printf("★★★★★★★★★★★★★★★★★★★★★★★★★\n");
		printf("★                                              ★\n");
		printf("★             알고리즘 프로그램                ★\n");
		printf("★                                              ★\n");
		printf("★★★★★★★★★★★★★★★★★★★★★★★★★\n\n"); 
		printf("=================================================\n");
		printf("l              0.메인 메뉴                       l\n");
		printf("=================================================\n");
		printf("l              1.정렬 알고리즘                   l\n");
		printf("=================================================\n");
		printf("l              2.그래프 탐색 알고리즘            l\n");
		printf("=================================================\n");
		printf("l              3.최소신장트리 알고리즘           l\n");
		printf("=================================================\n");
		printf("l              4.상태공간트리 알고리즘           l\n");
		printf("=================================================\n");
		printf("l              5.종료                            l\n");
		printf("=================================================\n\n");

		printf("▶ 번호를 입력하세요: ");
		scanf_s("%d", &number);
		
		system("cls");		// 콘솔 화면 초기화 + #include<windows.h> 헤더파일 포함

		if (number == 1) {
			

			printf("\n▷ 정렬 알고리즘을 선택하셨습니다.\n");
			printf("=================================================\n");
			printf("l                  0.메인 메뉴                  l\n");
			printf("=================================================\n");
			printf("l                  1.선택정렬                   l\n");
			printf("=================================================\n");
			printf("l                  2.버블정렬                   l\n");
			printf("=================================================\n");
			printf("l                  3.삽입정렬                   l\n");
			printf("=================================================\n");
			printf("l                  4.병합정렬                   l\n");
			printf("=================================================\n");
			printf("l                  5.퀵정렬                     l\n");
			printf("=================================================\n");
			printf("l                  6.힙정렬                     l\n");
			printf("=================================================\n");
			printf("l                  7.기수정렬                   l\n");
			printf("=================================================\n");
			printf("l                  8.계수정렬                   l\n");
			printf("=================================================\n");
			printf("l                  9.종료                       l\n");
			printf("=================================================\n\n");

			printf("▶ 번호를 입력하세요: ");
			scanf_s("%d", &sortnumber);
			if (0 <= sortnumber && sortnumber < 9) {
				system("cls");		// 콘솔 화면 초기화
				switch (sortnumber) {
					case 0:
						printf("▷ 메인 메뉴로 돌아가길 선택하셨습니다.\n");
						break;
					case 1:
						printf("\n▷ 선택정렬을 선택하였습니다.\n");
						setting();
						selectSort(arr, size);
						termination();
						break;
					case 2:
						printf("\n▷ 버블정렬을 선택하였습니다.\n");
						setting();
						BubbleSort(arr, size);
						termination();
						break;
					case 3:
						printf("\n▷ 삽입정렬을 선택하였습니다.\n");
						setting();
						InsertionSort(arr, size);
						termination();
						break;
					case 4:
						printf("\n▷ 병합정렬을 선택하였습니다.\n");
						setting();
						mergeSort(arr, 0, size - 1);
						termination();
						break;
					case 5:
						printf("\n▷ 퀵정렬을 선택하였습니다.\n");
						setting();
						quickSort(arr, 0, size - 1, size);
						termination();
						break;
					case 6:
						printf("\n▷ 힙정렬을 선택하였습니다.\n");
						setting();

						printf("정렬 전>> \n");
						printArray(arr, size);
						HeapSort(arr, size);
						printf("정렬 후 >>\n");
						printArray(arr, size);
						termination();
						break;
					case 7:
						printf("\n▷ 기수정렬을 선택하였습니다.\n");
						setting();

						printf("정렬 전>>\n");
						printArray(arr, size);
						radixSort(arr, size);
						printf("정렬 후 >>\n");
						printArray(arr, size);
						termination();
						break;
					case 8:
						printf("\n▷ 계수정렬을 선택하였습니다.\n");
						setting();

						printf("정렬 전>>\n");
						for (int i = 0; i < size; i++) {
							printf("%d ", arr[i]);
						}
						printf("\n");
						countingSort(arr, size);
						termination();
						break;
					case 9:
						printf("▷ 프로그램을 종료합니다.\n");
						return 0;
				}
			}

			else {
				printf("\n");
				printf("▷ 프로그램을 종료합니다.\n");
				return 0;
			}

		}

		else if (number == 2) {
			int graphnum;

			printf("\n▷ 그래프 알고리즘을 선택하셨습니다.\n");

			printf("=================================================\n");
			printf("l                  0.메인 메뉴                  l\n");
			printf("=================================================\n");
			printf("l                  1.깊이우선탐색               l\n");
			printf("=================================================\n");
			printf("l                  2.너비우선탐색               l\n");
			printf("=================================================\n");
			printf("l                  3.종료                       l\n");
			printf("=================================================\n\n");

			printf("▶ 번호를 입력하세요: ");
			scanf_s("%d", &graphnum);
			system("cls");		// 콘솔 화면 초기화
			//printf("%d번을 선택하셨습니다.\n", graphnum);

			if (graphnum == 1) {
				printf("▷ 깊이 우선 탐색을 시작합니다.\n");
				int i;
				graphType* G9;
				G9 = (graphType*)malloc(sizeof(graphType));
				createGraph(G9);

				for (i = 0; i < 7; i++)
					insertVertex(G9, i);
				insertEdge(G9, 0, 2);
				insertEdge(G9, 0, 1);
				insertEdge(G9, 1, 4);
				insertEdge(G9, 1, 3);
				insertEdge(G9, 1, 0);
				insertEdge(G9, 2, 4);
				insertEdge(G9, 2, 0);
				insertEdge(G9, 3, 6);
				insertEdge(G9, 3, 1);
				insertEdge(G9, 4, 6);
				insertEdge(G9, 4, 2);
				insertEdge(G9, 4, 1);
				insertEdge(G9, 5, 6);
				insertEdge(G9, 6, 5);
				insertEdge(G9, 6, 4);
				insertEdge(G9, 6, 3);
				printf("\n G9의 인접 리스트: ");
				print_adjList(G9);

				printf("\n\n/////////////\n\n깊이 우선 탐색>> ");
				DFS_adjList(G9, 0);

				getchar();
				termination();
			}

			else if (graphnum == 2) {
				printf("▷ 너비 우선 탐색을 시작합니다.\n");
				int i;
				graphType* G9;
				G9 = (graphType*)malloc(sizeof(graphType));
				createGraph(G9);

				for (i = 0; i < 7; i++)
					insertVertex(G9, i);
				insertEdge(G9, 0, 2);
				insertEdge(G9, 0, 1);
				insertEdge(G9, 1, 4);
				insertEdge(G9, 1, 3);
				insertEdge(G9, 1, 0);
				insertEdge(G9, 2, 4);
				insertEdge(G9, 2, 0);
				insertEdge(G9, 3, 6);
				insertEdge(G9, 3, 1);
				insertEdge(G9, 4, 6);
				insertEdge(G9, 4, 2);
				insertEdge(G9, 4, 1);
				insertEdge(G9, 5, 6);
				insertEdge(G9, 6, 5);
				insertEdge(G9, 6, 4);
				insertEdge(G9, 6, 3);
				printf("\n G9의 인접 리스트: ");
				print_adjList(G9);

				printf("\n\n/////////////\n\n깊이 우선 탐색>> ");
				BFS_adjList(G9, 0);

				getchar(); 
				termination();

			}
			else if (graphnum == 0) {
				printf("▷ 메인 메뉴로 돌아가길 선택하셨습니다.\n");
				continue;
			}

			else if (graphnum == 3) {
				printf("▷ 프로그램을 종료합니다.\n");
				return 0;
			}
			else {
				printf("▷ 프로그램을 종료합니다.");
				return 0;
			}
		}

		else if (number == 3)
		{
			int mstnum;
			printf("\n▷ 최소신장트리 알고리즘을 선택하셨습니다.\n");
			
			printf("================================================\n");
			printf("l               0.메인 메뉴                    l\n");
			printf("================================================\n");
			printf("l               1.프림 알고리즘                l\n");
			printf("===============================================\n");
			printf("l               2.크루스칼 알고리즘            l\n");
			printf("================================================\n");
			printf("l               3.종료                         l\n");
			printf("================================================\n\n");

			printf("▶ 번호를 입력하세요: ");

			scanf_s("%d", &mstnum);
			system("cls");		// 콘솔 화면 초기화
			//printf("%d번을 선택하셨습니다.\n", mstnum);

			if (mstnum == 1) {
				printf("\n▷ 프림 알고리즘을 시작합니다.\n");
				// 인접 행렬로 표현된 예시 그래프
				int G[MAX_VERTICES][MAX_VERTICES] = {
					{0, 2, 0, 6, 0},
					{2, 0, 3, 8, 5},
					{0, 3, 0, 0, 7},
					{6, 8, 0, 0, 9},
					{0, 5, 7, 9, 0}
				};

				int V = 5; // 정점 수
				int startVertex = 0; // 시작 정점
				printf("\n인접 행렬로 표현된 예시 그래프\n\t{ 0, 2, 0, 6, 0 },\n\t{ 2, 0, 3, 8, 5 },\n\t{ 0, 3, 0, 0, 7 },\n\t{ 6, 8, 0, 0, 9 },\n\t{ 0, 5, 7, 9, 0 }\n");
				printf("\nV(정점 수) : 5\nr(시작 정점) : 0\n");
				printf("\n");
				prim(G, V, startVertex);

				termination();
			}
			else if (mstnum == 2) {
				printf("\n▷ 크루스칼 알고리즘을 시작합니다.\n");
				int V = 4; // 정점 수
				int E = 5; // 간선 수
				Edge edges[MAX_EDGES] = {
					{0, 1, 10},
					{0, 2, 6},
					{0, 3, 5},
					{1, 3, 15},
					{2, 3, 4}
				};
				printf("\n인접 행렬로 표현된 예시 그래프\n\t{ 0, 1, 10 },\n\t{ 0, 2, 6 },\n\t{ 0, 3, 5 },\n\t{ 1, 3, 15 },\n\t{ 2, 3, 4 }\n");
				printf("\nV(정점 수) : 4\nE(간선 수) : 5\n");
				kruskal(edges, V, E);

				termination();
			}
			else if (mstnum == 0) {
				printf("▷ 메인 메뉴로 돌아가기를 선택하셨습니다.\n");
			}
			else if (mstnum == 3) {
				printf("▷ 프로그램을 종료합니다.\n");
				return 0;
			}

			else {
				printf("▷ 프로그램을 종료합니다.\n");
			}
		}
		else if (number == 0) {
			system("cls");
			printf("▷ 메인 메뉴로 돌아가기를 선택하셨습니다.\n");
			continue;
		}
		else if (number == 4) {
			int mstnum;
			printf("\n▷ 상태공간트리 알고리즘을 선택하셨습니다.\n");

			printf("================================================\n");
			printf("l               0.메인 메뉴                    l\n");
			printf("================================================\n");
			printf("l               1.A* 알고리즘                  l\n");
			printf("===============================================\n");
			printf("l               2.백트래킹                     l\n");
			printf("================================================\n");
			printf("l               3.종료                         l\n");
			printf("================================================\n\n");

			printf("▶ 번호를 입력하세요: ");

			scanf_s("%d", &mstnum);
			system("cls");		// 콘솔 화면 초기화

			if (mstnum == 0) {
				printf("▷ 메인 메뉴로 돌아가기를 선택하셨습니다.\n");
			}
			else if (mstnum == 1) {
				printf("▷ A* 알고리즘을 시작합니다.\n");
				int grid[ROW][COL] = {
				   {1, 1, 1, 1, 1},
				   {0, 0, 1, 0, 1},
				   {1, 1, 1, 1, 1},
				   {1, 0, 0, 0, 1},
				   {1, 1, 1, 1, 1}
				};

				printf("***  경로  ***\n");
				printf("■ ■ ■ ■ ■\n");
				printf("□ □ ■ □ ■\n");
				printf("■ ■ ■ ■ ■\n");
				printf("■ □ □ □ ■\n");
				printf("■ ■ ■ ■ ■\n\n");

				Point src = { 0, 0 };
				Point dest = { 4, 4 };

				printf("시작점: (%d, %d)\n", src.row, src.col);
				printf("목적지: (%d, %d)\n\n", dest.row, dest.col);

				aStarSearch(grid, src, dest);

				termination();
			}
			else if (mstnum == 2) {
				printf("▷ 백트래킹 알고리즘을 시작합니다.\n\n");

				int n, r;

				printf("▷ 원소의 개수를 입력하세요 (n): ");
				scanf_s("%d", &n);

				printf("▷ 조합의 크기를 입력하세요 (r): ");
				scanf_s("%d", &r);

				printf("▷ 배열의 원소를 입력하세요 (구분자 Enter):\n");
				for (int i = 0; i < n; ++i) {
					scanf_s("%d", &inputArray[i]);
				}

				combination(n, r, 0);

				termination();
			}
			else if (mstnum == 3){
				printf("▷ 프로그램을 종료합니다.");
				return 0;
			}
		}
		else if (number == 5) {
			printf("▷ 프로그램을 종료합니다.\n");
			return 0;
		}

		else {
			printf("▷ 번호를 잘못 입력하였습니다. 다시 시도해주세요.\n");
			continue;
		}
	}
}