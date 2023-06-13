#!/bin/bash

for x in *.xyz; do
    cat top.com > ${x%}.com
    cat ${x} >> ${x%}.com
    cat bottom.com >> ${x%}.com
    sed -i -e '8,9d' ${x%}.com
    echo "" >> ${x%}.com
done
