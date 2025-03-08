setErrorExitMsg() {
    set -E
    trap 'echo "ERROR: [$(basename "$0")] Failed at line $LINENO" >&2' ERR
}

cdSilent(){
    args="$@"
    cd "$args" > /dev/null
}