import os
import sys
import platform
import shutil
import time

try:
    import pymake
except:
    msg = 'Error. Pymake package is not available.\n'
    msg += 'Try installing using the following command:\n'
    msg += ' pip install https://github.com/modflowpy/pymake/zipball/master'
    raise Exception(msg)

try:
    import flopy
except:
    msg = 'Error. FloPy package is not available.\n'
    msg += 'Try installing using the following command:\n'
    msg += ' pip install flopy'
    raise Exception(msg)

import targets

sfmt = '{:25s} - {}'

class Simulation(object):

    def __init__(self, name):
        delFiles = True
        for idx, arg in enumerate(sys.argv):
            if arg.lower() == '--keep':
                if len(sys.argv) > idx + 1:
                    delFiles = False
                    break

        msg = sfmt.format('Initializing test', name)
        print(msg)
        self.name = name
        self.simpath = None
        self.inpt = None
        self.outp = None

        sysinfo = platform.system()
        self.delFiles = delFiles
        if sysinfo.lower() == 'windows':
            self.delFiles = False
        self.success = False
        return

    def __repr__(self):
        return self.name

    def setup(self, src, dst):
        msg = sfmt.format('Setup test', self.name)
        print(msg)
        self.originpath = src
        self.simpath = dst
        # write message
        print('running pymake.setup_mf6 from ' +
              '{}'.format(os.path.abspath(os.getcwd())))
        try:
            self.inpt, self.outp = pymake.setup_mf6(src=src, dst=dst)
            print('waiting...')
            time.sleep(0.5)
            success = True
        except:
            success = False
            print('source:      {}'.format(src))
            print('destination: {}'.format(dst))
        assert success, 'did not run pymake.setup_mf6'

        # Copy comparison simulations if available
        if success:
            action = pymake.setup_mf6_comparison(src, dst,
                                                 remove_existing=self.delFiles)

            self.action = action
        return

    def run(self):
        """
        Run the model and assert if the model terminated successfully
        """
        msg = sfmt.format('Run test', self.name)
        print(msg)

        # Set nam as namefile name without path
        nam = 'mfsim.nam'

        # run mf6 models
        exe = os.path.abspath(targets.target_dict[targets.program])
        msg = sfmt.format('using executable', exe)
        print(msg)
        try:
            success, buff = flopy.run_model(exe, nam, model_ws=self.simpath,
                                            silent=False, report=True)
            msg = sfmt.format('MODFLOW 6 run', self.name)
            if success:
                print(msg)
            else:
                print(msg)
        except:
            msg = sfmt.format('MODFLOW 6 run', self.name)
            print(msg)
            success = False

        assert success

        self.nam_cmp = None
        if success:
            if self.action is not None:
                if self.action.lower() == 'compare':
                    msg = sfmt.format('Comparison files', self.name)
                    print(msg)
                else:
                    cpth = os.path.join(self.simpath, self.action)
                    key = self.action.lower().replace('.cmp', '')
                    exe = os.path.abspath(targets.target_dict[key])
                    npth = pymake.get_namefiles(cpth)[0]
                    nam = os.path.basename(npth)
                    self.nam_cmp = nam
                    try:
                        success_cmp, buff = flopy.run_model(exe, nam,
                                                            model_ws=cpth,
                                                            silent=False,
                                                            report=True)
                        msg = sfmt.format('Comparison run',
                                          self.name + '/' + key)
                        if success:
                            print(msg)
                        else:
                            print(msg)
                    except:
                        success_cmp = False
                        msg = sfmt.format('Comparison run',
                                          self.name + '/' + key)
                        print(msg)

                    assert success_cmp

        return

    def compare(self):
        """
        Compare the model results

        """
        msg = sfmt.format('Comparison test', self.name)
        print(msg)
        
        success_tst = False
        if self.action is not None:
            cpth = os.path.join(self.simpath, self.action)
            files_cmp = None
            if self.action.lower() == 'compare':
                files_cmp = []
                files = os.listdir(cpth)
                for file in files:
                    files_cmp.append(
                        os.path.abspath(os.path.join(cpth, file)))

            files1 = []
            files2 = []
            exfiles = []
            ipos = 0
            for file1 in self.outp:
                ext = os.path.splitext(file1)[1]
                if ext.lower() == '.hds':
                    pth = os.path.join(self.simpath, file1)
                    files1.append(pth)
                    if files_cmp is None:
                        files2 = None
                    else:
                        pth = os.path.join(cpth, file1 + '.cmp')
                        files2.append(pth)
                        txt = sfmt.format('Comparison file {}'.format(ipos+1),
                                          os.path.basename(pth))
                        print(txt)

                    # look for an exclusion file
                    pth = os.path.join(self.simpath, file1 + '.ex')
                    if os.path.isfile(pth):
                        exfiles.append(pth)
                    else:
                        exfiles.append(None)

                    # increment ipos
                    ipos += 1

            if self.nam_cmp is None:
                pth = None
            else:
                pth = os.path.join(cpth, self.nam_cmp)

            for ipos in range(len(files1)):
                file1 = files1[ipos]
                outfile = os.path.splitext(os.path.basename(file1))[0]
                outfile = os.path.join(self.simpath, outfile + '.hds.cmp.out')
                if files2 is None:
                    file2 = None
                else:
                    file2 = files2[ipos]

                # set exfile
                exfile = None
                if file2 is None:
                    if len(exfiles) > 0:
                        exfile = exfiles[ipos]
                        if exfile is not None:
                            txt = sfmt.format('Exclusion file {}'.format(ipos + 1),
                                              os.path.basename(exfile))
                            print(txt)

                # make comparison
                success_tst = pymake.compare_heads(None, pth,
                                                   precision='double',
                                                   outfile=outfile,
                                                   files1=file1,
                                                   files2=file2,
                                                   verbose=True,
                                                   exfile=exfile)
                msg = sfmt.format('Head comparison {}'.format(ipos+1), 
                                  self.name)
                if success_tst:
                    print(msg)
                else:
                    print(msg)

                assert success_tst, msg

        self.success = success_tst
        return

    def teardown(self):
        """
        Remove the example folder

        """
        if self.success:
            if self.delFiles:
                msg = sfmt.format('Teardown test', self.name)
                print(msg)
                try:
                    shutil.rmtree(self.simpath)
                    success = True
                except:
                    print('Could not remove test ' + self.name)
                    success = False
                assert success
            else:
                print('Retaining test files')
        return
