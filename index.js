const express = require('express');  // require express
const { spawn } = require('child_process');  // child process module allows to execute any shell commands or scripts within nodejs
const upload = require('express-fileupload');  // module allows file upload

const app = express();  // create the app
var port = 20236
var host = '0.0.0.0'
app.listen(port, host, () => console.log(`Listening at port ${port}... Go to http://${host}:${port}/`));  //listen at a port (port number 3000)
app.use(express.static('public'));

app.use(upload())

app.get('/', (req, res) => {
    res.sendFile(__dirname + 'index.html')
})

app.post('/', (req, res) => {
    if (req.files) {
        console.log(req.files)
        var file = req.files.file
        var filename = file.name
        console.log(filename)

        file.mv('./uploads/'+filename, function (err) {
            if (err) {
                res.send(err)
            } else {
                res.send("File Uploaded")
            }
        })
    }
})

// const ls = spawn('python3.9', ['temp/parse_mmcif_to_domains.py', 'temp/AF-P43238-F1-model_v2.cif']);

// ls.stdout.on('data', (data) => {
//   console.log(`stdout: ${data}`);
// });

// ls.stderr.on('data', (data) => {
//   console.error(`stderr: ${data}`);
// });

// ls.on('close', (code) => {
//   console.log(`child process exited with code ${code}`);
// }); 