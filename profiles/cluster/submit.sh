#!/usr/bin/env bash

select_account() {
  if [[ -n "${SNAKEMAKE_SLURM_ACCOUNT:-}" ]]; then
    echo "${SNAKEMAKE_SLURM_ACCOUNT}"
    return 0
  fi

  local acct=""
  acct="$(sacctmgr --noheader --parsable2 show user "$USER" format=DefaultAccount 2>/dev/null | head -n1 | tr -d '[:space:]' || true)"
  if [[ -n "${acct}" && "${acct}" != "Unknown" ]]; then
    echo "$acct"
    return 0
  fi

  # Otherwise pick the first account returned
  acct="$(sacctmgr --noheader --parsable2 show user "$USER" format=Account 2>/dev/null | head -n1 | tr -d '[:space:]' || true)"
  if [[ -n "${acct}" && "${acct}" != "Unknown" ]]; then
    echo "$acct"
    return 0
  fi

  echo ""
}

have_account_arg=0
for a in "$@"; do
  [[ "$a" == --account* || "$a" == -A* ]] && have_account_arg=1 && break
done

acct="$(select_account)"
if [[ $have_account_arg -eq 0 && -n "${acct}" ]]; then
  exec sbatch --account="$acct" "$@"
else
  exec sbatch "$@"
fi
