// #define SS_DSDPDATA_DUMP
#include "dsdphsd.h"

#include "env/copt_env.h"
#include "base/copt_prob.h"
#include "base/copt_qmatrix.h"

#include "tool/tool_assert.h"
#include "tool/tool_define.h"
#include "tool/tool_memory.h"
#include "tool/tool_index.h"
#include "tool/tool_math.h"
#include "tool/tool_sort.h"
#include "tool/tool_time.h"

#include "lp/cone_solver.h"

int main(int argc, char* argv[])
{
  ss_retcode retcode = RETCODE_OK;
  copt_env* env = NULL;
  copt_prob* prob = NULL;

  if (argc != 2)
  {
    printf("\n%s <sdpafile>\n", argv[0]);
    return 0;
  }

  SS_CALL(ss_Env_Init(NULL, NULL, &env));
  SS_CALL(ss_Prob_Create(env, &prob));

  SS_CALL(ss_Prob_ReadSDPA(prob, argv[1]));
  SS_CALL(ss_Prob_SolveHDSDP(prob));

exit_cleanup:
  if (retcode != RETCODE_OK)
  {
    printf("Error when solving with HDSDP.\n");
  }

  ss_Prob_Delete(&prob);
  ss_Env_Free(&env);
  return 0;
}
