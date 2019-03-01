#1/bin/bash
#set -x
MESH_DIR=${PWD}/data
LOG_DIR=${PWD}/log/negSumDiag
TEST_FILE=test.cpp

if [ -d ${MESH_DIR} ]; then
	cd ${MESH_DIR}
else
	echo "Error: Directory ${MESH_DIR} does not exist!"
	exit 1
fi

OWNER=`find ${MESH_DIR} -maxdepth 1 -name "owner*" | sort`
NEIGHBOR=`find ${MESH_DIR} -maxdepth 1 -name "neighbour*" | sort`
OWNER_ARR=(${OWNER})
NEIGHBOR_ARR=(${NEIGHBOR})
length=${#NEIGHBOR_ARR[@]}
echo "length: ${length}"
cd ../
DIM_ARR=(1 3 6 9)
dims=${#DIM_ARR[@]}
for ((i=0;i<length;i++))
do
	OWNER_NAME=`echo ${OWNER_ARR[i]} | sed 's/.*\/\([^\/]\)/\1/g' `
	NEIGHBOR_NAME=`echo ${NEIGHBOR_ARR[i]} | sed 's/.*\/\([^\/]\)/\1/g' `
	NONZERONUM=`cat ${NEIGHBOR_ARR[i]} | grep nInternalFaces | sed 's/.*nInternalFaces:\ \([0-9]*\).*/\1/g' `
	MESH_SIZE=`echo ${OWNER_NAME} | sed 's/owner_\(.*\)/\1/g' `
	echo "Start test on case ${MESH_SIZE}"
#	echo $OWNER_NAME
#	echo ${NEIGHBOR_NAME}
	echo ${NONZERONUM}
	`sed -i "s/owner\[\]\ =\ \".*\"/owner\[\]\ =\ \"data\/${OWNER_NAME}\"/g" ${TEST_FILE} `
	`sed -i "s/neighbor\[\]\ =\ \".*\"/neighbor\[\]\ =\ \"data\/${NEIGHBOR_NAME}\"/g" ${TEST_FILE} `
	`sed -i "s/NONZERONUM\ \([0-9]*\)/NONZERONUM\ ${NONZERONUM}/g" ${TEST_FILE}`
	for ((j=0;j<dims;j++))
	do
		`sed -i "s/DIMENSION\ \([0-9]*\)/DIMENSION\ ${DIM_ARR[j]}/g" ${TEST_FILE}`
		`sh mkrun > ${LOG_DIR}/log_${MESH_SIZE}_${DIM_ARR[j]}`
		SLAVE=`cat ${LOG_DIR}/log_${MESH_SIZE}_${DIM_ARR[j]} | grep Slave | sed 's/Slave Core Time:\ \([0-9]\)/\1/g' `
		MASTER=`cat ${LOG_DIR}/log_${MESH_SIZE}_${DIM_ARR[j]} | grep CPU | sed 's/CPU Processor Time:\ \([0-9]\)/\1/g' `
		CORRECT=`cat ${LOG_DIR}/log_${MESH_SIZE}_${DIM_ARR[j]} | grep correct`
		echo "Mesh size: ${MESH_SIZE} Dimension: ${DIM_ARR[j]} Slave: ${SLAVE} Master: ${MASTER} Result: ${CORRECT}"
	done
#	tmp=`cat ${TEST_FILE} | grep "char owner" | sed "s/\".*\"/\"${OWNER_NAME}\"/g" `
done
#echo ${OWNER}
