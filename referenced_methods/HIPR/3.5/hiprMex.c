/* Maximum flow - highest lavel push-relabel algorithm */
/* COPYRIGHT C 1995, 2000 by IG Systems, Inc., igsys@eclipse.net */
#include "mex.h"
#include "types.h"  /* type definitions */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

// #include "types.h"  /* type definitions */
// #include "parser2.c" /* parser */
/*
#include "timer.c"           timing routine 
*/

/*
#define GLOB_UPDT_FREQ 0.5
*/
#define GLOB_UPDT_FREQ 0
#define ALPHA 6
#define BETA 12
#define MAXLONG 100000

#define MAXLINE       10000	/* max line length in the input file */
#define ARC_FIELDS      3	/* no of fields in arc line  */
#define NODE_FIELDS     2	/* no of fields in node line  */
#define P_FIELDS        3       /* no of fields in problem line */
#define PROBLEM_TYPE "max"      /* name of problem type*/

#define WHITE 0
#define GREY 1
#define BLACK 2

/* global variables */

long n;          /* number of nodes */
long m;          /* number of arcs */
long nm;         /* n + ALPHA * m */
long nMin;       /* smallest node id */
node *nodes;     /* array of nodes */
arc *arcs;       /* array of arcs */
bucket *buckets; /* array of buckets */
cType *cap;      /* array of capacities */
node *source;    /* source node pointer */
node *sink;      /* sink node pointer */
//node   **queue;              /* queue for BFS */
//node   **qHead, **qTail, **qLast;     /* queue pointers */
long dMax;                /* maximum label */
long aMax;                /* maximum actie node label */
long aMin;                /* minimum active node label */
double flow;              /* flow value */
long pushCnt = 0;         /* number of pushes */
long relabelCnt = 0;      /* number of relabels */
long updateCnt = 0;       /* number of updates */
long checkCnt = 0;        /* number of check node and discharge */
long gapCnt = 0;          /* number of gaps */
long gNodeCnt = 0;        /* number of nodes after gap */
float t, t2;              /* for saving times */
node *sentinelNode;       /* end of the node list marker */
arc *stopA;               /* used in forAllArcs */
long workSinceUpdate = 0; /* the number of arc scans since last update */
float globUpdtFreq;       /* global update frequency */

/* macros */

#define forAllNodes(i) for (i = nodes; i != sentinelNode; i++)
#define forAllArcs(i, a) for (a = i->first, stopA = (i + 1)->first; a != stopA; a++)

#define nNode(i) ((i)-nodes + nMin)
#define nArc(a) ((a == NULL) ? -1 : (a)-arcs)

#define min(a, b) (((a) < (b)) ? a : b)

long i_dist;

#define aAdd(l, i)             \
  {                            \
    i->bNext = l->firstActive; \
    l->firstActive = i;        \
    i_dist = i->d;             \
    if (i_dist < aMin)         \
      aMin = i_dist;           \
    if (i_dist > aMax)         \
      aMax = i_dist;           \
    if (dMax < aMax)           \
      dMax = aMax;             \
  }

/* i must be the first element */
#define aRemove(l, i)          \
  {                            \
    l->firstActive = i->bNext; \
  }

node *i_next, *i_prev;
#define iAdd(l, i)             \
  {                            \
    i_next = l->firstInactive; \
    i->bNext = i_next;         \
    i->bPrev = sentinelNode;   \
    i_next->bPrev = i;         \
    l->firstInactive = i;      \
  }

#define iDelete(l, i)               \
  {                                 \
    i_next = i->bNext;              \
    if (l->firstInactive == i)      \
    {                               \
      l->firstInactive = i_next;    \
      i_next->bPrev = sentinelNode; \
    }                               \
    else                            \
    {                               \
      i_prev = i->bPrev;            \
      i_prev->bNext = i_next;       \
      i_next->bPrev = i_prev;       \
    }                               \
  }

/* allocate datastructures, initialize related variables */

int allocDS()

{

  nm = ALPHA * n + m; /* global relabelling parameter   */
  /*
  queue = (node**) calloc ( n, sizeof (node*) );
  if ( queue == NULL ) return ( 1 );
  qLast = queue + n - 1;
  qInit();
  */
  buckets = (bucket *)calloc(n + 2, sizeof(bucket));
  if (buckets == NULL)
    return (1);

  sentinelNode = nodes + n;
  sentinelNode->first = arcs + 2 * m;

  return (0);

} /* end of allocate */

void init()

{
  node *i; /* current node */
  int overflowDetected;
  bucket *l;
  arc *a;
#ifdef EXCESS_TYPE_LONG
  double testExcess;
#endif
#ifndef OLD_INIT
  unsigned long delta;
#endif

  // initialize excesses

  forAllNodes(i)
  {
    i->excess = 0;
    i->current = i->first;
    forAllArcs(i, a)
        a->resCap = cap[a - arcs];
  }

  for (l = buckets; l <= buckets + n - 1; l++)
  {
    l->firstActive = sentinelNode;
    l->firstInactive = sentinelNode;
  }

  overflowDetected = 0;
#ifdef EXCESS_TYPE_LONG
  testExcess = 0;
  forAllArcs(source, a)
  {
    if (a->head != source)
    {
      testExcess += a->resCap;
    }
  }
  if (testExcess > MAXLONG)
  {
    printf("c WARNING: excess overflow. See README for details.\nc\n");
    overflowDetected = 1;
  }
#endif
#ifdef OLD_INIT
  source->excess = MAXLONG;
#else
  if (overflowDetected)
  {
    source->excess = MAXLONG;
  }
  else
  {
    source->excess = 0;
    forAllArcs(source, a)
    {
      if (a->head != source)
      {
        pushCnt++;         /* push + 1*/
        delta = a->resCap; /* flow */
        a->resCap -= delta;
        (a->rev)->resCap += delta;
        a->head->excess += delta; /* update */
      }
    }
  }

  /*  setup labels and buckets */
  l = buckets + 1;

  aMax = 0;
  aMin = n;

  forAllNodes(i)
  {
    if (i == sink)
    {
      i->d = 0;
      iAdd(buckets, i);
      continue;
    }
    if ((i == source) && (!overflowDetected))
    {
      i->d = n;
    }
    else
      i->d = 1;
    if (i->excess > 0)
    {
      /* put into active list */
      aAdd(l, i);
    }
    else
    { /* i -> excess == 0 */
      /* put into inactive list */
      if (i->d < n)
        iAdd(l, i);
    }
  }
  dMax = 1;
#endif

  //  dMax = n-1;
  //  flow = 0.0;

} /* end of init */

void checkMax()

{
  bucket *l;

  for (l = buckets + dMax + 1; l < buckets + n; l++)
  {
    assert(l->firstActive == sentinelNode);
    assert(l->firstInactive == sentinelNode);
  }
}

/* global update via backward breadth first search from the sink */

void globalUpdate()

{

  node *i, *j;    /* node pointers */
  arc *a;         /* current arc pointers  */
  bucket *l, *jL; /* bucket */
  long curDist, jD;
  long state;

  updateCnt++;

  /* initialization */

  forAllNodes(i)
      i->d = n;
  sink->d = 0;

  for (l = buckets; l <= buckets + dMax; l++)
  {
    l->firstActive = sentinelNode;
    l->firstInactive = sentinelNode;
  }

  dMax = aMax = 0;
  aMin = n;

  /* breadth first search */

  // add sink to bucket zero

  iAdd(buckets, sink);
  for (curDist = 0; 1; curDist++)
  {

    state = 0;
    l = buckets + curDist;
    jD = curDist + 1;
    jL = l + 1;
    /*
    jL -> firstActive   = sentinelNode;
    jL -> firstInactive  = sentinelNode;
    */

    if ((l->firstActive == sentinelNode) &&
        (l->firstInactive == sentinelNode))
      break;

    while (1)
    {

      switch (state)
      {
      case 0:
        i = l->firstInactive;
        state = 1;
        break;
      case 1:
        i = i->bNext;
        break;
      case 2:
        i = l->firstActive;
        state = 3;
        break;
      case 3:
        i = i->bNext;
        break;
      default:
        assert(0);
        break;
      }

      if (i == sentinelNode)
      {
        if (state == 1)
        {
          state = 2;
          continue;
        }
        else
        {
          assert(state == 3);
          break;
        }
      }

      /* scanning arcs incident to node i */
      forAllArcs(i, a)
      {
        if (a->rev->resCap > 0)
        {
          j = a->head;
          if (j->d == n)
          {
            j->d = jD;
            j->current = j->first;
            if (jD > dMax)
              dMax = jD;

            if (j->excess > 0)
            {
              /* put into active list */
              aAdd(jL, j);
            }
            else
            {
              /* put into inactive list */
              iAdd(jL, j);
            }
          }
        }
      } /* node i is scanned */
    }
  }

} /* end of global update */

/* second stage -- preflow to flow */
void stageTwo()
/*
   do dsf in the reverse flow graph from nodes with excess
   cancel cycles if found
   return excess flow in topological order
*/

/*
   i->d is used for dfs labels 
   i->bNext is used for topological order list
   buckets[i-nodes]->firstActive is used for DSF tree
*/

{
  node *i, *j, *tos, *bos, *restart, *r;
  arc *a;
  cType delta;

  /* deal with self-loops */
  forAllNodes(i)
  {
    forAllArcs(i, a) if (a->head == i)
    {
      a->resCap = cap[a - arcs];
    }
  }

  /* initialize */
  tos = bos = NULL;
  forAllNodes(i)
  {
    i->d = WHITE;
    //    buckets[i-nodes].firstActive = NULL;
    buckets[i - nodes].firstActive = sentinelNode;
    i->current = i->first;
  }

  /* eliminate flow cycles, topologicaly order vertices */
  forAllNodes(i) if ((i->d == WHITE) && (i->excess > 0) &&
                     (i != source) && (i != sink))
  {
    r = i;
    r->d = GREY;
    do
    {
      for (; i->current != (i + 1)->first; i->current++)
      {
        a = i->current;
        if ((cap[a - arcs] == 0) && (a->resCap > 0))
        {
          j = a->head;
          if (j->d == WHITE)
          {
            /* start scanning j */
            j->d = GREY;
            buckets[j - nodes].firstActive = i;
            i = j;
            break;
          }
          else if (j->d == GREY)
          {
            /* find minimum flow on the cycle */
            delta = a->resCap;
            while (1)
            {
              delta = min(delta, j->current->resCap);
              if (j == i)
                break;
              else
                j = j->current->head;
            }

            /* remove delta flow units */
            j = i;
            while (1)
            {
              a = j->current;
              a->resCap -= delta;
              a->rev->resCap += delta;
              j = a->head;
              if (j == i)
                break;
            }

            /* backup DFS to the first saturated arc */
            restart = i;
            for (j = i->current->head; j != i; j = a->head)
            {
              a = j->current;
              if ((j->d == WHITE) || (a->resCap == 0))
              {
                j->current->head->d = WHITE;
                if (j->d != WHITE)
                  restart = j;
              }
            }

            if (restart != i)
            {
              i = restart;
              i->current++;
              break;
            }
          }
        }
      }

      if (i->current == (i + 1)->first)
      {
        /* scan of i complete */
        i->d = BLACK;
        if (i != source)
        {
          if (bos == NULL)
          {
            bos = i;
            tos = i;
          }
          else
          {
            i->bNext = tos;
            tos = i;
          }
        }

        if (i != r)
        {
          i = buckets[i - nodes].firstActive;
          i->current++;
        }
        else
          break;
      }
    } while (1);
  }

  /* return excesses */
  /* note that sink is not on the stack */
  if (bos != NULL)
  {
    for (i = tos; i != bos; i = i->bNext)
    {
      a = i->first;
      while (i->excess > 0)
      {
        if ((cap[a - arcs] == 0) && (a->resCap > 0))
        {
          if (a->resCap < i->excess)
            delta = a->resCap;
          else
            delta = i->excess;
          a->resCap -= delta;
          a->rev->resCap += delta;
          i->excess -= delta;
          a->head->excess += delta;
        }
        a++;
      }
    }
    /* now do the bottom */
    i = bos;
    a = i->first;
    while (i->excess > 0)
    {
      if ((cap[a - arcs] == 0) && (a->resCap > 0))
      {
        if (a->resCap < i->excess)
          delta = a->resCap;
        else
          delta = i->excess;
        a->resCap -= delta;
        a->rev->resCap += delta;
        i->excess -= delta;
        a->head->excess += delta;
      }
      a++;
    }
  }
}

/* gap relabeling */

int gap(emptyB)
    bucket *emptyB;

{

  bucket *l;
  node *i;
  long r; /* index of the bucket before l  */
  int cc; /* cc = 1 if no nodes with positive excess before
		      the gap */

  gapCnt++;
  r = (emptyB - buckets) - 1;

  /* set labels of nodes beyond the gap to "infinity" */
  for (l = emptyB + 1; l <= buckets + dMax; l++)
  {
    /* this does nothing for high level selection 
    for (i = l -> firstActive; i != sentinelNode; i = i -> bNext) {
      i -> d = n;
      gNodeCnt++;
    }
    l -> firstActive = sentinelNode;
    */

    for (i = l->firstInactive; i != sentinelNode; i = i->bNext)
    {
      i->d = n;
      gNodeCnt++;
    }

    l->firstInactive = sentinelNode;
  }

  cc = (aMin > r) ? 1 : 0;

  dMax = r;
  aMax = r;

  return (cc);
}

/*--- relabelling node i */

long relabel(i)

    node *i; /* node to relabel */

{

  node *j;
  long minD; /* minimum d of a node reachable from i */
  arc *minA; /* an arc which leads to the node with minimal d */
  arc *a;

  assert(i->excess > 0);

  relabelCnt++;
  workSinceUpdate += BETA;

  i->d = minD = n;
  minA = NULL;

  /* find the minimum */
  forAllArcs(i, a)
  {
    workSinceUpdate++;
    if (a->resCap > 0)
    {
      j = a->head;
      if (j->d < minD)
      {
        minD = j->d;
        minA = a;
      }
    }
  }

  minD++;

  if (minD < n)
  {

    i->d = minD;
    i->current = minA;

    if (dMax < minD)
      dMax = minD;

  } /* end of minD < n */

  return (minD);

} /* end of relabel */

/* discharge: push flow out of i until i becomes inactive */

void discharge(i)

    node *i;

{

  node *j;    /* sucsessor of i */
  long jD;    /* d of the next bucket */
  bucket *lj; /* j's bucket */
  bucket *l;  /* i's bucket */
  arc *a;     /* current arc (i,j) */
  cType delta;
  arc *stopA;

  assert(i->excess > 0);
  assert(i != sink);

  checkCnt++;

  do
  {

    jD = i->d - 1;
    l = buckets + i->d;

    /* scanning arcs outgoing from  i  */
    for (a = i->current, stopA = (i + 1)->first; a != stopA; a++)
    {
      if (a->resCap > 0)
      {
        j = a->head;

        if (j->d == jD)
        {
          pushCnt++;
          if (a->resCap < i->excess)
            delta = a->resCap;
          else
            delta = i->excess;
          a->resCap -= delta;
          a->rev->resCap += delta;

          if (j != sink)
          {

            lj = buckets + jD;

            if (j->excess == 0)
            {
              /* remove j from inactive list */
              iDelete(lj, j);
              /* add j to active list */
              aAdd(lj, j);
            }
          }

          j->excess += delta;
          i->excess -= delta;

          if (i->excess == 0)
            break;

        } /* j belongs to the next bucket */
      }   /* a  is not saturated */
    }     /* end of scanning arcs from  i */

    if (a == stopA)
    {
      /* i must be relabeled */
      checkCnt++;
      relabel(i);

      if (i->d == n)
        break;
      if ((l->firstActive == sentinelNode) &&
          (l->firstInactive == sentinelNode))
        gap(l);

      if (i->d == n)
        break;
    }
    else
    {
      /* i no longer active */
      i->current = a;
      /* put i on inactive list */
      iAdd(l, i);
      break;
    }
  } while (1);
}

// go from higher to lower buckets, push flow
void wave()
{

  node *i;
  bucket *l;

  for (l = buckets + aMax; l > buckets; l--)
  {
    for (i = l->firstActive; i != sentinelNode; i = l->firstActive)
    {
      aRemove(l, i);

      assert(i->excess > 0);
      discharge(i);
    }
  }
}

/* first stage  -- maximum preflow*/

void stageOne()

{

  node *i;
  bucket *l; /* current bucket */

#if defined(INIT_UPDATE) || defined(OLD_INIT) || defined(WAVE_INIT)
  globalUpdate();
#endif

  workSinceUpdate = 0;

#ifdef WAVE_INIT
  wave();
#endif

  /* main loop */
  while (aMax >= aMin)
  {
    l = buckets + aMax; /* highest label */
    i = l->firstActive; /* First in, last out? */

    if (i == sentinelNode)
      aMax--;
    else
    {
      aRemove(l, i);

      assert(i->excess > 0);
      discharge(i);

      if (aMax < aMin)
        break;

      /* is it time for global update? */
      if (workSinceUpdate * globUpdtFreq > nm)
      {
        globalUpdate();
        workSinceUpdate = 0;
      }
    }

  } /* end of the main loop */

  flow = sink->excess;
}

int parse2(long *n_ad, long *m_ad, node **nodes_ad, arc **arcs_ad, long **cap_ad, node **source_ad, node **sink_ad, long *node_min_ad ,const mxArray* prhs[]){



  long    n,                      /* internal number of nodes */
          node_min=0,             /* minimal no of node  */
          node_max=0,             /* maximal no of nodes */
        *arc_first = NULL,         /* internal array for holding
                                      - node degree
                                      - position of the first outgoing arc */
        *arc_tail = NULL,          /* internal array: tails of the arcs */
          source=0,               /* no of the source */
          sink=0,                 /* no of the sink */
          /* temporary variables carrying no of nodes */
          head, tail, i;

  long    m,                      /* internal number of arcs */
          /* temporary variables carrying no of arcs */
          last, arc_num, arc_new_num;

  node    *nodes=NULL,            /* pointer to the node structure */
          *head_p,
          *ndp;

  arc     *arcs=NULL,             /* pointer to the arc structure */
          *arc_current=NULL,
          *arc_new,
          *arc_tmp;

  long    *acap=NULL,             /* array of capasities */
          cap;                    /* capasity of the current arc */

  long    no_lines=0,             /* no of current input line */
          no_plines=0,            /* no of problem-lines */
          no_nslines=0,           /* no of node-source-lines */
          no_nklines=0,           /* no of node-source-lines */
          no_alines=0,            /* no of arc-lines */
          pos_current=0;          /* 2*no_alines */

  char    in_line[MAXLINE],       /* for reading input line */
          pr_type[3],             /* for reading type of the problem */
          nd;                     /* source (s) or sink (t) */

  int     k,                      /* temporary */
          err_no;                 /* no of detected error */

  /* -------------- error numbers & error messages ---------------- */
  #define EN1   0
  #define EN2   1
  #define EN3   2
  #define EN4   3
  #define EN6   4
  #define EN10  5
  #define EN7   6
  #define EN8   7
  #define EN9   8
  #define EN11  9
  #define EN12 10
  #define EN13 11
  #define EN14 12
  #define EN16 13
  #define EN15 14
  #define EN17 15
  #define EN18 16
  #define EN21 17
  #define EN19 18
  #define EN20 19
  #define EN22 20

  static char *err_message[] = 
    { 
  /* 0*/    "more than one problem line.",
  /* 1*/    "wrong number of parameters in the problem line.",
  /* 2*/    "it is not a Max Flow problem line.",
  /* 3*/    "bad value of a parameter in the problem line.",
  /* 4*/    "can't obtain enough memory to solve this problem.",
  /* 5*/    "more than one line with the problem name.",
  /* 6*/    "can't read problem name.",
  /* 7*/    "problem description must be before node description.",
  /* 8*/    "this parser doesn't support multiply sources and sinks.",
  /* 9*/    "wrong number of parameters in the node line.",
  /*10*/    "wrong value of parameters in the node line.",
  /*11*/    " ",
  /*12*/    "source and sink descriptions must be before arc descriptions.",
  /*13*/    "too many arcs in the input.",
  /*14*/    "wrong number of parameters in the arc line.",
  /*15*/    "wrong value of parameters in the arc line.",
  /*16*/    "unknown line type in the input.",
  /*17*/    "reading error.",
  /*18*/    "not enough arcs in the input.",
  /*19*/    "source or sink doesn't have incident arcs.",
  /*20*/    "can't read anything from the input file."
    };
  /* --------------------------------------------------------------- */

  n = mxGetM(prhs[0])+2;
  m = mxGetM(prhs[1])*2;
  double *NN = mxGetPr(prhs[2]);
  double *ss = mxGetPr(prhs[0]);
  double *ee = mxGetPr(prhs[1]);
  int N = NN[0];
  int mEE = mxGetM(prhs[1]);
  int mSS = mxGetM(prhs[0]);
  m = m + 2*N;

  /* alloc */
  nodes    = (node*) calloc ( n+2, sizeof(node) );
  arcs     = (arc*)  calloc ( 2*m+1, sizeof(arc) );
  arc_tail = (long*) calloc ( 2*m,   sizeof(long) ); 
  arc_first= (long*) calloc ( n+2, sizeof(long) );
  acap     = (long*) calloc ( 2*m, sizeof(long) );
  arc_current = arcs;
  /* source and sink */
  source = n-1;
  sink = n;
  node_max = 0;
  node_min = n;

  /* edge */
  for(int line = 0;line<mEE;line++){
      /* node 1 -> node 2 */
      //////////////////////
      tail = ee[line];
      head = ee[line + mEE];
      cap = ee[line + mEE*2];
      arc_first[tail + 1] ++; 
      arc_first[head + 1] ++;

      /* storing information about the arc */
      arc_tail[pos_current]        = tail;
      arc_tail[pos_current+1]      = head;
      arc_current       -> head    = nodes + head;
      arc_current       -> resCap    = cap;
      arc_current       -> rev  = arc_current + 1;
      ( arc_current + 1 ) -> head    = nodes + tail;
      ( arc_current + 1 ) -> resCap    = 0;
      ( arc_current + 1 ) -> rev  = arc_current;

      /* searching minimumu and maximum node */
      if ( head < node_min ) node_min = head;
      if ( tail < node_min ) node_min = tail;
      if ( head > node_max ) node_max = head;
      if ( tail > node_max ) node_max = tail;

      no_alines   ++;
      arc_current += 2;
      pos_current += 2;

      /* node 2 -> node 1*/
      //////////////////////
      long tmp = tail;
      tail = head;
      head = tmp;
      cap = ee[line + mEE*3];
      arc_first[tail + 1] ++; 
      arc_first[head + 1] ++;

      /* storing information about the arc */
      arc_tail[pos_current]        = tail;
      arc_tail[pos_current+1]      = head;
      arc_current       -> head    = nodes + head;
      arc_current       -> resCap    = cap;
      arc_current       -> rev  = arc_current + 1;
      ( arc_current + 1 ) -> head    = nodes + tail;
      ( arc_current + 1 ) -> resCap    = 0;
      ( arc_current + 1 ) -> rev  = arc_current;

      /* searching minimumu and maximum node */
      if ( head < node_min ) node_min = head;
      if ( tail < node_min ) node_min = tail;
      if ( head > node_max ) node_max = head;
      if ( tail > node_max ) node_max = tail;

      no_alines   ++;
      arc_current += 2;
      pos_current += 2;
  }

  for(int line = 0;line<mSS;line++){
    if(ss[line]>0){
        /////////////////////////////
        tail = source;
        head = line + 1;
        cap = ss[line];

        arc_first[tail + 1] ++; 
        arc_first[head + 1] ++;

        /* storing information about the arc */
        arc_tail[pos_current]        = tail;
        arc_tail[pos_current+1]      = head;
        arc_current       -> head    = nodes + head;
        arc_current       -> resCap    = cap;
        arc_current       -> rev  = arc_current + 1;
        ( arc_current + 1 ) -> head    = nodes + tail;
        ( arc_current + 1 ) -> resCap    = 0;
        ( arc_current + 1 ) -> rev  = arc_current;

        /* searching minimumu and maximum node */
        if ( head < node_min ) node_min = head;
        if ( tail < node_min ) node_min = tail;
        if ( head > node_max ) node_max = head;
        if ( tail > node_max ) node_max = tail;

        no_alines   ++;
        arc_current += 2;
        pos_current += 2;
    }
    if(ss[line+mSS]>0){
        ////////////////////////////
        tail = line+1;
        head = sink;
        cap = ss[tail+mSS-1];
        arc_first[tail + 1] ++; 
        arc_first[head + 1] ++;

        /* storing information about the arc */
        arc_tail[pos_current]        = tail;
        arc_tail[pos_current+1]      = head;
        arc_current       -> head    = nodes + head;
        arc_current       -> resCap    = cap;
        arc_current       -> rev  = arc_current + 1;
        ( arc_current + 1 ) -> head    = nodes + tail;
        ( arc_current + 1 ) -> resCap    = 0;
        ( arc_current + 1 ) -> rev  = arc_current;

        /* searching minimumu and maximum node */
        if ( head < node_min ) node_min = head;
        if ( tail < node_min ) node_min = tail;
        if ( head > node_max ) node_max = head;
        if ( tail > node_max ) node_max = tail;

        no_alines   ++;
        arc_current += 2;
        pos_current += 2;
    }
  }

  /********** ordering arcs - linear time algorithm ***********/

  /* first arc from the first node */
  ( nodes + node_min ) -> first = arcs;

  /* before below loop arc_first[i+1] is the number of arcs outgoing from i;
    after this loop arc_first[i] is the position of the first 
    outgoing from node i arcs after they would be ordered;
    this value is transformed to pointer and written to node.first[i]
    */
  
  for ( i = node_min + 1; i <= node_max + 1; i ++ ) 
    {
      arc_first[i]          += arc_first[i-1];          /* summation of arc till now */
      ( nodes + i ) -> first = arcs + arc_first[i];     /* node->first, the address of first arc */
    }


  for ( i = node_min; i < node_max; i ++ ) /* scanning all the nodes  
                                              exept the last*/
                                              /* Just like sort the arc, put all arcs in the position they should be */
    {

      last = ( ( nodes + i + 1 ) -> first ) - arcs;       /* The last arc for node i */
                              /* arcs outgoing from i must be cited    
                                from position arc_first[i] to the position
                                equal to initial value of arc_first[i+1]-1  */

      for ( arc_num = arc_first[i]; arc_num < last; arc_num ++ )
        { tail = arc_tail[arc_num];     /* tail of current arc */

    while ( tail != i )
            /* the arc no  arc_num  is not in place because arc cited here
              must go out from i;
              we'll put it to its place and continue this process
              until an arc in this position would go out from i */

      { arc_new_num  = arc_first[tail];
        arc_current  = arcs + arc_num;
        arc_new      = arcs + arc_new_num;
        
        /* arc_current must be cited in the position arc_new    
          swapping these arcs:                                 */

        head_p               = arc_new -> head;
        arc_new -> head      = arc_current -> head;
        arc_current -> head  = head_p;

        cap                 = arc_new -> resCap;
        arc_new -> resCap     = arc_current -> resCap;
        arc_current -> resCap = cap;

        if ( arc_new != arc_current -> rev )
          {
            arc_tmp                = arc_new -> rev;
            arc_new  -> rev     = arc_current -> rev;
            arc_current -> rev  = arc_tmp;

                  ( arc_current -> rev ) -> rev = arc_current;
      ( arc_new     -> rev ) -> rev = arc_new;
          }

        arc_tail[arc_num] = arc_tail[arc_new_num];
        arc_tail[arc_new_num] = tail;

        /* we increase arc_first[tail]  */
        arc_first[tail] ++ ;

              tail = arc_tail[arc_num];
      }
        }
      /* all arcs outgoing from  i  are in place */
    }       

  /* -----------------------  arcs are ordered  ------------------------- */

  /*----------- constructing lists ---------------*/
  /* rearrange the arc list to each node */

    for ( ndp = nodes + node_min; ndp <= nodes + node_max;  ndp ++ )
        ndp -> first = (arc*) NULL;
  /* descending order, to save the first pointer of arc for node ndp */
    for ( arc_current = arcs + (2*m-1); arc_current >= arcs; arc_current -- )
      {
        arc_num = arc_current - arcs;
        tail = arc_tail [arc_num];
        ndp = nodes + tail;
        /* avg
        arc_current -> next = ndp -> first;
        */
        ndp -> first = arc_current;
      }


  /* ----------- assigning output values ------------*/
  *m_ad = m;
  *n_ad = node_max - node_min + 1;
  *source_ad = nodes + source;
  *sink_ad   = nodes + sink;
  *node_min_ad = node_min;
  *nodes_ad = nodes + node_min;
  *arcs_ad = arcs;
  *cap_ad = acap;

  for ( arc_current = arcs, arc_num = 0; 
        arc_num < 2*m;
        arc_current ++, arc_num ++
      )
      acap [ arc_num ] = arc_current -> resCap; 

  if ( source < node_min || source > node_max )
    /* bad value of the source */
    { err_no = EN20; goto error; }
    
  if ( (*source_ad) -> first == (arc*) NULL ||
      (*sink_ad  ) -> first == (arc*) NULL    ) 
    /* no arc goes out of the source */
    { err_no = EN20; goto error; }

  /* free internal memory */
  free ( arc_first ); free ( arc_tail );

  /* Thanks God! all is done */
  return (0);

  /* ---------------------------------- */
  error:  /* error found reading input */

  printf ( "\nline %ld of input - %s\n", 
          no_lines, err_message[err_no] );

}
/* --------------------   end of parser  -------------------*/


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  int cc;
  globUpdtFreq = GLOB_UPDT_FREQ;
  checkCnt = 0;
  parse2(&n, &m, &nodes, &arcs, &cap, &source, &sink, &nMin, prhs);

  cc = allocDS();
  if (cc)
  {
    fprintf(stderr, "Allocation error\n");
  }

  /*  t = timer();
  t2 = t; */

  init();
  stageOne();

  /* t2 = timer() - t2; */

  // printf("c flow:       %12.01f\n", flow);

  #ifndef CUT_ONLY
    stageTwo();

    /* t = timer() - t; */

    // printf("c time:        %10.2f\n", t);

  #endif

    // printf("c cut tm:      %10.2f\n", t2);

  #ifdef CHECK_SOLUTION

  #endif

  free ( nodes-1 ); 
  free ( arcs ); 
  free ( cap ); 
  free( buckets );

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *maxFlow = mxGetPr(plhs[0]);
  maxFlow[0] = flow;

  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double* pushNum = mxGetPr(plhs[1]);
  pushNum[0] = checkCnt;

  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double* updateNum = mxGetPr(plhs[2]);
  updateNum[0] = updateCnt;
  // plhs[0] = mxCreateDoubleScalar(1);
};
