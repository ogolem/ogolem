pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
                echo 'Building  ogolem and running unit tests'
                sh 'gradle build -i'
            }
        }
    }
}
